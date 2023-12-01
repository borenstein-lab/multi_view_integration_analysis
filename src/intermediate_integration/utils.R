###################################
# Utilities
###################################

require(BiocParallel)
require(cowplot)
require(logger)

is_valid_r_name <- function(string) {
  grepl("^((([[:alpha:]]|[.][._[:alpha:]])[._[:alnum:]]*)|[.])$", string)
}

remove_samples_from_input <- function(diablo_input, held_out_samples) {
  new_diablo_input <- list()
  new_diablo_input$X <- lapply(diablo_input$X, function(m) {m[-held_out_samples,]})
  new_diablo_input$Y <- diablo_input$Y[-held_out_samples]
  new_diablo_input$sample_ids <- diablo_input$sample_ids[-held_out_samples]
  return(new_diablo_input)
}

# Gets a list of p values and returns significance stars
get_signif_mark <- function(ps) {
  sapply(ps, function(p) {
    if (p < 0.0001) return("***")
    if (p < 0.001) return("**")
    if (p < 0.05) return("*")
    if (p < 0.1) return(".")
    return("")
  })
}

# Function for organizing data given in one wide table into an input suitable 
#  for running the 'DIABLO' function.
# The single table 'proc_data' is expected to have one sample_id column, one 
#  label column, and features divided into 'views', each view represented by a 
#  prefix to the feature name in the form of a single letter and two underscores 
#  (e.g. 'T__' for taxonomy, 'P__' for pathways).
# Note: MintTea was not tested for more than 4 views.
organize_data_for_diablo <- function(
    proc_data, 
    study_group_column, 
    sample_id_column, 
    view_prefixes,
    min_features_per_view = 5) {
  
  # Verify that sample_id and study_group columns are present
  if (! study_group_column %in% colnames(proc_data)) stop("Could not find study_group_column in data table") 
  if (! sample_id_column %in% colnames(proc_data)) stop("Could not find sample_id_column in data table")
  ds_label <- proc_data[[study_group_column]] 
  ds_sample_ids <- proc_data[[sample_id_column]]
  
  # Verify all feature names are valid names
  feat_names <- colnames(proc_data)[! colnames(proc_data) %in% c(study_group_column, sample_id_column)]
  if (any(! sapply(feat_names, is_valid_r_name)))
    stop('Found invalid feature names. See: https://www.w3schools.com/r/r_variables_name.asp')
  
  # Verify all features indeed start with one of the provided prefixes
  is_prefix_valid <- function(feat_name, prefixes) { any(startsWith(feat_name, prefix = prefixes)) }
  valid_prefix_check <- sapply(feat_names, is_prefix_valid, prefixes = paste0(view_prefixes, '__'))
  if (any(!valid_prefix_check)) 
    stop('Found features that do not start with any of the given view-prefixes (followed by two underscores)')
  
  # Organize features by groups of views. Check sufficient number of features. 
  ds_X_list <- list()
  for (prfx in view_prefixes) {
    tmp <- proc_data %>% dplyr::select(starts_with(paste0(prfx, "__"))) %>% as.matrix()
    if (ncol(tmp) < min_features_per_view) 
      stop("Not enough features in view: '", prfx, 
           "'. Expecting a minimum of ", min_features_per_view)
    ds_X_list[[prfx]] <- tmp
  }
  
  return(list(
    X = ds_X_list,
    Y = ds_label,
    sample_ids = ds_sample_ids
  ))
}

# Given a table with features from different types, return pairwise correlations between features from different omics only
get_all_inter_omic_corrs <- function(comp_data) {
  all_corrs_pears <- cor(comp_data) %>% 
    abs() %>%
    data.frame() %>% 
    tibble::rownames_to_column('feature1') %>% 
    tidyr::pivot_longer(cols = -feature1, names_to = "feature2", values_to = "pearson_corr") 
  
  all_corrs_spear <- cor(comp_data, method = 'spearman') %>% 
    abs() %>%
    data.frame() %>% 
    tibble::rownames_to_column('feature1') %>% 
    tidyr::pivot_longer(cols = -feature1, names_to = "feature2", values_to = "spearman_corr")

  all_corrs <- left_join(all_corrs_pears, 
                         all_corrs_spear, 
                         by = c('feature1', 'feature2')) %>%
    filter(substr(feature1, 1, 1) != substr(feature2, 1, 1)) %>%
    filter(feature1 > feature2)
  
  return(all_corrs)
}



get_auc <- function(pred_col, data, label_col = "true_label", auto_direction = FALSE) {
  if (n_distinct(data[[label_col]]) == 1) return(NA)
  f <- as.formula(paste(label_col, "~", pred_col))
  if (auto_direction) return(as.numeric(auc(roc(formula = f, data = data, quiet = TRUE))))
  return(as.numeric(auc(roc(formula = f, data = data, levels = c("healthy", "disease"), direction = "<"))))
}

get_cv_auc <- function(diablo_input, rf_cv_results = NULL, folds, list_keepX, final_ncomp, diablo_dist, diablo_design) {
  require(pROC)
  require(rsample)
  
  n_samples <- length(diablo_input$sample_ids)
  
  # Collect AUC stats
  my_aucs <- data.frame()
  
  for (i in 1:nrow(folds)) {
    fold_id <- paste(folds$splits[[i]]$id, collapse = " ")
    cat('.')
    
    # Partition the data into train and test ("held-out")
    # held_out_samples_i <- sample(1:n_samples, round(n_samples/10))
    # train_samples <- setdiff(1:n_samples, held_out_samples_i)
    
    train_samples <- folds$splits[[i]]$in_id
    held_out_samples_i <- setdiff(1:n_samples, train_samples)
    new_diablo_input <- remove_samples_from_input(diablo_input, held_out_samples_i)
    held_out_data <- remove_samples_from_input(diablo_input, train_samples)
    
    # Train DIABLO model
    tmp_diablo_out <- block.splsda(
      new_diablo_input$X, 
      new_diablo_input$Y, 
      ncomp = final_ncomp, 
      keepX = list_keepX, 
      max.iter = 300,
      design = diablo_design
    )
    
    # Make prediction on held out samples
    predictions <- predict(tmp_diablo_out, held_out_data$X, dist = diablo_dist)
    
    # Calculate auroc
    tmp <- data.frame(preds = predictions$AveragedPredict[,1,], true_label = held_out_data$Y)
    #tmp2 <- bind_rows(tmp2, tmp %>% mutate(fold_id=fold_id))
    my_aucs_v <- lapply(1:final_ncomp, function(i) get_auc(paste0("preds.dim", i), tmp))
    names(my_aucs_v) <- paste0("auc_comp", 1:final_ncomp)
    my_aucs <- bind_rows(
      my_aucs, 
      data.frame(my_aucs_v) %>%
        mutate(fold_id = fold_id) %>%
        mutate(h_ratio_train = round(unname(table(new_diablo_input$Y)['healthy']) / length(new_diablo_input$Y), 3)) %>%
        mutate(h_ratio_test = round(unname(table(held_out_data$Y)['healthy']) / length(held_out_data$Y), 3))
    )
  }
  cat('\n')
  
  # Add information from ML pipeline results (matched by fold)
  if (!is.null(rf_cv_results))
    my_aucs <- my_aucs %>%
      left_join(
        rf_cv_results %>% 
          select(fold_id, out_of_fold_test_auc, mean_out_of_fold_test_auc, h_ratio_train, h_ratio_test) %>% 
          rename(rf_fold_auc = out_of_fold_test_auc) %>% 
          rename(mean_rf_auc = mean_out_of_fold_test_auc) %>% 
          rename(h_ratio_train_SANITY = h_ratio_train, h_ratio_test_SANITY = h_ratio_test),
        by = "fold_id"
      )
  
  return(my_aucs)
}

get_cor_btwn_variates <- function(diablo_out, block1, block2, comp_id) {
  comp_txts <- paste0('comp', comp_id)
  if (length(comp_txts) == 1) 
    comp_txts <- c(comp_txts)
  unname(sapply(
    comp_txts, 
    function(comp_txt) {
      round(cor(diablo_out$variates[[block1]][,comp_txt], 
                diablo_out$variates[[block2]][,comp_txt]), 3) 
  }))
}

# get_ttest_variate_vs_label <- function(diablo_out, block1, comp_id) {
#   comp_txts <- paste0('comp', comp_id)
#   if (length(comp_txts) == 1) comp_txts <- c(comp_txts)
#   unname(sapply(comp_txts, function(comp_txt) {
#     variate <- diablo_out$variates[[block1]][,comp_txt]
#     t.test(variate[diablo_out$Y == 'healthy'], variate[diablo_out$Y != 'healthy'])$p.value 
#   }))
# }
# 
# Re-organizes the loadings from a diablo model object into a single compact table
organize_diablo_loadings <- function(diablo_out) {
  loadings_per_block <- lapply(
    1:(length(diablo_out$loadings)-1), # Ignore the last table (dummy loadings for Y)
    function(i) {
      df <- diablo_out$loadings[[i]]
      comp_cols <- colnames(df)
      df <- df %>%
        as.data.frame() %>%
        tibble::rownames_to_column("feature") %>%
        mutate(feature_set = names(diablo_out$loadings)[i]) %>%
        tidyr::pivot_longer(cols = all_of(comp_cols),
                            names_to = "component",
                            values_to = "loading",
                            names_prefix = "comp") %>%
        mutate(component = as.numeric(component)) %>%
        filter(loading != 0)
    })
  return(bind_rows(loadings_per_block))
}

# organize_cv_feature_stability <- function(diablo_out, diablo_cv) {
#   stab_per_repeat <- bind_rows(
#     # Iterate over the number of CV repeats
#     lapply(
#       1:length(diablo_cv$features$stable), 
#       function(k) 
#         bind_rows(
#           # Iterate over the omics
#           lapply(
#             1:length(diablo_out$X), 
#             function(j) 
#               bind_rows(
#                 # Iterate over the components
#                 lapply(
#                   1:unname(diablo_out$ncomp[1]), 
#                   function(i) 
#                     data.frame(diablo_cv$features$stable[[k]][[j]][[i]]) %>% 
#                       mutate(Var1 = as.character(Var1)) %>% 
#                       rename(feature = Var1, feature_cv_stability = Freq) %>% 
#                       mutate(component = i))) %>% 
#                       mutate(feature_set = names(diablo_out$X)[j]))) %>% 
#                       mutate(nrep = k)
#     )
#   )
#   
#   # Summarize over repeats
#   n_reps <- max(stab_per_repeat$nrep)
#   mean_stab <- stab_per_repeat %>%
#     group_by(component, feature_set, feature) %>%
#     summarise(
#       # Note that I divide by the number of repeats to consider repeats where no 
#       #  fold selected the features will also be taken into account as 0 stability
#       mean_cv_stability = sum(feature_cv_stability) / n_reps, 
#       .groups = "drop")
#   return(mean_stab)
# }
# 
# organize_diablo_loadings_with_stability <- function(diablo_out, diablo_cv) {
#   diablo_loadings <- organize_diablo_loadings(diablo_out)
#   stability <- organize_cv_feature_stability(diablo_out, diablo_cv)
#   return(diablo_loadings %>% 
#            left_join(stability, 
#                      by = c('component', 'feature_set', 'feature')) %>%
#            mutate(mean_cv_stability = 
#                     tidyr::replace_na(mean_cv_stability, 0)))
# }
# 
# # A method to test the significance of each canonical component, inspired by 
# #  the permutation approach in Witten D.M. et al, 2009.
# # See also 'MultiCCA.permute' function in 'PMA' package. 
# # Roughly, this approach checks whether our components detect inter-omic
# #  correlations that are higher than those expected by random.
# # To expand this method to the supervised case as in DIABLO, we also check
# #  the model's error rate for each permutation.
# 
# test_significance_perms <- function(
#     diablo_out, 
#     diablo_dist, 
#     n_perms = 99, 
#     permute_label = FALSE, 
#     permute_omics = TRUE) {
#   
#   n_ds <- length(diablo_out$X)
#   n_comp <- unname(diablo_out$ncomp['Y'])
#   cors_per_comp <- data.frame()
#   err_per_comp <- data.frame()
#   
#   # Run permutations (perm 0 = no permutation at all)
#   set.seed(27)
#   for(i in 0:n_perms) { 
#     if (i%%5==0) cat('.')
#     
#     # Permute samples / label only
#     if (i > 0) {
#       if (permute_omics) {
#         tmp_ds_list <- lapply(
#           diablo_out$X, 
#           function(ds) { 
#             ds_perm <- ds[sample(nrow(ds)),] 
#             rownames(ds_perm) <- as.character(1:nrow(ds))
#             return(ds_perm)
#             }
#           )
#       }
#       if (permute_label) {
#         tmp_label <- sample(diablo_out$Y)
#       }
#     } else {
#       tmp_ds_list <- diablo_out$X
#       tmp_label <- diablo_out$Y
#     }
#     
#     # Train model with permuted data
#     tmp_diablo_out <- block.splsda(
#       tmp_ds_list, 
#       tmp_label, 
#       ncomp = n_comp, 
#       keepX = diablo_out$keepX, 
#       max.iter = 200,
#       design = diablo_out$design
#     )
#     
#     # Evaluate predictive power of the model (BER)
#     tmp_diablo_cv <- perf(
#       tmp_diablo_out, 
#       validation = 'Mfold', 
#       folds = 10, 
#       nrepeat = 10,  
#       progressBar = FALSE,
#       dist = diablo_dist,
#       auc = FALSE
#       )
#     
#     # Record BER for this permutation
#     err_per_comp <- bind_rows(
#       err_per_comp,
#       data.frame(
#         perm = i,
#         comp = 1:n_comp,
#         #auc = unname(sapply(tmp_diablo_cv$auc, function(x) {x[,'AUC']})),
#         BER = unname(tmp_diablo_cv$WeightedPredict.error.rate['Overall.BER',])
#       )
#     )
#     # ggplot(err_per_comp, aes(x = perm, y = BER, color = (perm==0))) + geom_point(size = 3) + facet_wrap(~ comp) + theme_classic()
#     
#     # Iterate over pairs of feature sets to compute between-omic correlations
#     for (ds1 in 1:(n_ds-1)) {
#       for (ds2 in (ds1+1):n_ds) {
#         # Get and record correlations
#         cors_per_comp <- bind_rows(
#           cors_per_comp,
#           data.frame(
#             perm = i,
#             comp = 1:n_comp,
#             ds_pair = paste(ds1, ds2, sep='_'),
#             cor = get_cor_btwn_variates(tmp_diablo_out, ds1, ds2, 1:n_comp)
#           )
#         )
#       }
#     } # Completed loop over pairs of omics
#   } # Completed permutations loop
#   cat('\n')
#   
#   # Compute the sum of correlations per iteration and per component ID
#   sum_cor_per_comp <- cors_per_comp %>%
#     group_by(perm, comp) %>%
#     summarise(sum_cor = sum(cor), .groups = "drop")
#   true_sum_cors <- sum_cor_per_comp %>% filter(perm == 0)
#   sum_cor_per_comp <- sum_cor_per_comp %>% filter(perm > 0)
#   
#   # Get p-value per component, reflecting how likely are the *correlation values
#   #  between variates* in our model to be achieved by chance?
#   pvals_1 <- c()
#   for(i in 1:n_comp) {
#     n_perms_above_true <- sum(
#       sum_cor_per_comp %>% 
#           filter(comp == i) %>% 
#           pull(sum_cor) >= 
#         true_sum_cors %>% 
#           filter(comp == i) %>% 
#           pull(sum_cor)
#       )
#     pvals_1 <- c(
#       pvals_1, 
#       (n_perms_above_true + 1) / (n_perms + 1) # See North et al. doi: 10.1086/341527
#       )
#   }
#   
#   # Similarly, get p-value per component, reflecting how likely are the 
#   #  *error rates* of the overall model to be achieved by chance?
#   pvals_2 <- c()
#   for(i in 1:n_comp) {
#     n_perms_above_true <- sum(
#       err_per_comp %>% 
#         filter(comp == i) %>% 
#         pull(BER) >= 
#         (err_per_comp %>% 
#         filter(perm == 0) %>%
#         filter(comp == i) %>% 
#         pull(BER))
#     )
#     pvals_2 <- c(
#       pvals_2, 
#       (n_perms_above_true + 1) / (n_perms + 1) # See North et al. doi: 10.1086/341527
#     )
#   }
#   
#   # FDR and organize in table
#   res <- data.frame(
#     comp = 1:n_comp,
#     perm_BER_fdr = p.adjust(pvals_2, method = "fdr"),
#     perm_sum_cor_fdr = p.adjust(pvals_1, method = "fdr")
#   )
#   
#   return(res)
# }
###########################################
# Utilities for running/analyzing sGCCA
###########################################

get_loadings_smcca <- function(data_X, cca_out) {
  all_loadings <- data.frame()
  # Iterate over views
  for (k in 1:cca_out$K) {
    # Iterate over components
    for (comp in 1:ncol(cca_out$ws[[1]])) {
      # Collect non-zero loadings
      tmp <- data.frame(
        feature = colnames(data_X[[k]]),
        feature_set = names(data_X)[k],
        component = comp,
        ws = cca_out$ws[[k]][,comp]
      ) %>%
        filter(ws > 0 | ws < 0)
      all_loadings <- bind_rows(all_loadings, tmp)
    }
  }
  return(all_loadings)
}

###########################################
# Utilities for running/analyzing DIABLO
###########################################

# require(BiocParallel)
# require(cowplot)
# require(logger)
# 
# # Gets a list of p values and returns signifiance stars
# get_signif_mark <- function(ps) {
#   sapply(ps, function(p) {
#     if (p < 0.0001) return("***")
#     if (p < 0.001) return("**")
#     if (p < 0.05) return("*")
#     if (p < 0.01) return(".")
#     return("")
#   })
# }
# 
# # Function for organizing data for DIABLO for a given dataset - OLD
# organize_data_for_diablo <- function(d, proc_data) {
#   ds_label <- proc_data[[d]] %>% dplyr::pull(DiseaseState) 
#   ds_sample_ids <- proc_data[[d]] %>% dplyr::pull(sample_id__) 
#   ds_data_T <- proc_data[[d]] %>% dplyr::select(starts_with("T__")) %>% as.matrix()
#   ds_data_G <- proc_data[[d]] %>% dplyr::select(starts_with("G__")) %>% as.matrix() 
#   ds_data_P <- proc_data[[d]] %>% dplyr::select(starts_with("P__")) %>% as.matrix()
#   ds_data_M <- proc_data[[d]] %>% dplyr::select(starts_with("M__")) %>% as.matrix()
#   ds_list <- list("T" = ds_data_T, "G" = ds_data_G, "P" = ds_data_P, "M" = ds_data_M)  
#   
#   return(list(
#     X = ds_list,
#     Y = ds_label,
#     sample_ids = ds_sample_ids
#   ))
# }
# 
# # Function for organizing data for wrapped-DIABLO for a given dataset
# organize_data_for_diablo2 <- function(d, proc_data) {
#   ds_label <- proc_data %>% dplyr::pull(DiseaseState) 
#   ds_sample_ids <- proc_data %>% dplyr::pull(sample_id__) 
#   ds_data_T <- proc_data %>% dplyr::select(starts_with("T__")) %>% as.matrix()
#   ds_data_G <- proc_data %>% dplyr::select(starts_with("G__")) %>% as.matrix() 
#   ds_data_P <- proc_data %>% dplyr::select(starts_with("P__")) %>% as.matrix()
#   ds_data_M <- proc_data %>% dplyr::select(starts_with("M__")) %>% as.matrix()
#   no_metabs <- ifelse(ncol(ds_data_M)==0, TRUE, FALSE)
#   
#   ds_list <- list("T" = ds_data_T, "G" = ds_data_G, "P" = ds_data_P)  
#   if(!no_metabs) ds_list[['M']] <- ds_data_M
#     
#   return(list(
#     X = ds_list,
#     Y = ds_label,
#     sample_ids = ds_sample_ids
#   ))
# }
# 
# # Function for wrapped DIABLO (using specific pipeline parameters)
# run_wrapped_diablo <- function(
#     diablo_input,
#     diablo_design = 2/3,
#     n_repeats = 3,
#     n_folds = 10,
#     diablo_keepX = 10, 
#     diablo_ncomp = 5,
#     diablo_dist = 'centroids.dist',
#     seed = 27) {
#   
#   require(rsample)
#   require(mixOmics)
#   
#   if(! is.null(seed)) set.seed(seed) 
#   
#   # Set number of variables to take from each omic
#   list_keepX <- list()
#   for(x in names(diablo_input$X)) list_keepX[[x]] <- rep(diablo_keepX, diablo_ncomp)
#   
#   # Get number of samples
#   n_samples <- length(diablo_input$sample_ids)
#   
#   # Generate data subsets (=repeated cross-validation) 
#   folds <- vfold_cv(
#     data.frame(sample_num = 1:n_samples), 
#     v = n_folds, 
#     repeats = n_repeats
#   )
#   
#   # Place holders for AUC stats and component loadings
#   cv_loadings <- data.frame()
#   
#   for (i in 1:nrow(folds)) {
#     fold_id <- paste(folds$splits[[i]]$id, collapse = "_")
#     cat('.')
#     
#     # Partition the data into train and test ("held-out")
#     train_samples <- folds$splits[[i]]$in_id
#     held_out_samples_i <- setdiff(1:n_samples, train_samples)
#     new_diablo_input <- remove_samples_from_input(diablo_input, held_out_samples_i)
#     held_out_data <- remove_samples_from_input(diablo_input, train_samples)
#     
#     # Train DIABLO model
#     tmp_diablo_out <- block.splsda(
#       new_diablo_input$X, 
#       new_diablo_input$Y, 
#       ncomp = diablo_ncomp, 
#       keepX = list_keepX, 
#       max.iter = 300,
#       design = diablo_design
#     )
#     
#     # Make prediction on held out samples using DIABLO's approach 
#     # predictions <- predict(tmp_diablo_out, held_out_data$X, dist = diablo_dist)
#     
#     # Calculate auroc on held out data
#     # tmp <- data.frame(
#     #   preds = predictions$AveragedPredict[,1,], 
#     #   true_label = held_out_data$Y
#     #   )
#     # cv_aucs_v <- lapply(
#     #   1:diablo_ncomp, 
#     #   function(i) get_auc(paste0("preds.dim", i), tmp)
#     #   )
#     # names(cv_aucs_v) <- paste0("auc_comp", 1:diablo_ncomp)
#     # 
#     # cv_aucs <- bind_rows(
#     #   cv_aucs, 
#     #   data.frame(cv_aucs_v) %>%
#     #     mutate(fold_id = fold_id) 
#     # )
#     
#     # Save all non-zero loadings
#     tmp_loadings <- organize_diablo_loadings(tmp_diablo_out) %>%
#       # Add full component identifiers (i.e. component number + specific fold)
#       mutate(fold_id = fold_id) %>%
#       mutate(comp_fold_id = paste0(fold_id, '_Comp', component)) %>%
#       dplyr::select(-component)
#     cv_loadings <- bind_rows(cv_loadings, tmp_loadings)
#   }
#   cat('\n')
#   
#   # For each pair of features count the number of times they appeared together
#   # (Note - each pair currently appears twice: [A,B], [B,A])
#   feat_pairs <- 
#     # Get pairs of features that were found in the same component
#     cv_loadings %>%
#     full_join(cv_loadings, by = c("fold_id","comp_fold_id")) %>%
#     filter(feature.x != feature.y) %>%
#     # In cases where a pair of features come out together in more than one 
#     #  component in the same run, we only count it once
#     dplyr::select(feature.x, feature_set.x, feature.y, feature_set.y, fold_id) %>% 
#     distinct() %>%
#     group_by(feature.x, feature.y) %>%
#     summarise(N = n(), .groups = "drop") %>%
#     # Compute a distance between each pair, where 1 = never appeared together, 
#     #  i.e. maximal distance, 0 = always appear together, minimal distance
#     mutate(dist = ((n_repeats * n_folds) - N) / (n_repeats * n_folds)) 
#   
#   return(list(feat_pairs = feat_pairs))
# }
# 
# get_pairwise_cors_from_diablo_model <- function(diablo_out) {
#   blocks <- names(diablo_out$X)
#   n_blocks <- length(blocks)
#   n_comps <- unname(diablo_out$ncomp[1])
#   
#   cors <- data.frame()
#   
#   for (c in 1:n_comps) {
#     for (i in 1:(n_blocks-1)) {
#       for (j in (i+1):n_blocks) {
#         cur_cor <- cor.test(
#           diablo_out$variates[[blocks[i]]][,c],
#           diablo_out$variates[[blocks[j]]][,c]
#           )
#         cors <- bind_rows(
#           cors,
#           data.frame(
#             block1 = blocks[i],
#             block2 = blocks[j],
#             component = c,
#             cor = cur_cor$estimate,
#             cor.p.val = cur_cor$p.value
#           ))    
#       }
#     }
#   }
#   rownames(cors) <- NULL
#   return(cors)
# }

# Given a table with features from different types, return pairwise correlations between features from different omics only
get_all_inter_omic_corrs2 <- function(comp_data, control_for = NULL) {
  if (is.null(control_for)) {
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

  } else {
    control_for <- factor(control_for)
    all_corrs <- data.frame()
    for (ii in 2:ncol(comp_data)) {
      feature1 = colnames(comp_data)[ii]
      for (jj in 1:ii) {
        feature2 = colnames(comp_data)[jj]
        if (substr(feature1, 1, 1) != substr(feature2, 1, 1)) {
          tmp_cor <- tryCatch(pcor.test(comp_data[,ii], comp_data[,jj], as.numeric(control_for)),error=function(e) e, warning=function(w) w)
          tmp_sp_cor <- tryCatch(pcor.test(comp_data[,ii], comp_data[,jj], as.numeric(control_for), method = 'spearman'),error=function(e) e, warning=function(w) w)
          all_corrs <- bind_rows(
              all_corrs,
              data.frame(
                feature1 = feature1,
                feature2 = feature2,
                pearson_pcorr = ifelse(is.null(tmp_cor$estimate), 0, abs(tmp_cor$estimate)),
                spearman_pcorr = ifelse(is.null(tmp_sp_cor$estimate), 0, abs(tmp_sp_cor$estimate))
              )
            )
        }
      }
    }
    return(all_corrs)
  }
}

get_z_score_shared_links <- function(loadings_1, loadings_2, join_by = 'component', n = 100) {
  # Get common features and reduce both lists to common features only
  common_features <- intersect(loadings_1$feature, loadings_2$feature)
  if (length(common_features) <= 1) return(0)
  loadings_1 <- loadings_1 %>% filter(feature %in% common_features)
  loadings_2 <- loadings_2 %>% filter(feature %in% common_features)
  
  # Get all feature links in components #1
  feat_pairs_1 <- loadings_1 %>%
    full_join(loadings_1, by = c(join_by)) %>%
    filter(feature.x < feature.y) %>%
    mutate(feat_pair = paste(feature.x, feature.y)) %>%
    pull(feat_pair) 
  n1 <- n_distinct(feat_pairs_1)
  
  # Permutation test: How much more likely are features in components #2 to overlap with components #1, compared to shuffled components?
  perc_shared_links <- c()
  for (ii in 1:n) {
    tmp_2 <- loadings_2
    if (ii > 1) tmp_2[[join_by]] <- sample(tmp_2[[join_by]])
    feat_pairs_2 <- tmp_2 %>%
      full_join(tmp_2, by = c(join_by)) %>%
      filter(feature.x < feature.y) %>%
      mutate(feat_pair = paste(feature.x, feature.y)) %>%
      pull(feat_pair) 
    common_links <- intersect(feat_pairs_1, feat_pairs_2)
    perc_shared_links <- c(perc_shared_links, 100 * length(common_links) / n1)
  }
  if (sd(perc_shared_links[-1]) == 0) return(0)
  zscore = (perc_shared_links[1] - mean(perc_shared_links[-1])) / sd(perc_shared_links[-1])
  return(zscore)
}

# # Illustrates how the diablo design matrix impacts auc on the one side and 
# #  correlation between omic components on the other
# plot_correlation_discrimination_tradeoff <- function(
#     diablo_input, 
#     designs_to_try = seq(0,10,2) / 10, 
#     diablo_dist = "centroids.dist", 
#     diablo_prelim_ncomp = 5,
#     list_keepX = list_keepX,
#     n_repeats = 20,
#     BPPARAM = SerialParam()) 
# {
#   log_debug("Starting analysis of correlation-discrimination trade-off")
#   
#   # Collect AUC performance per design constant
#   auc_per_design <- lapply(designs_to_try, FUN = function(diablo_design) {
#     log_debug("Running DIABLO with design: {diablo_design}")
#     # Run DIABLO
#     diablo_out <- block.splsda(diablo_input$X, 
#                                diablo_input$Y, 
#                                ncomp = diablo_prelim_ncomp, 
#                                keepX = list_keepX, 
#                                design = diablo_design)
#     
#     cors <- get_pairwise_cors_from_diablo_model(diablo_out)
#     cors <- cors %>%
#       group_by(component) %>%
#       summarise(mean_pairwise_cor = mean(cor), sd_pairwise_cor = sd(cor)/sqrt(n()))
#     
#     # CV
#     diablo_cv <- perf(diablo_out, 
#                       validation = 'Mfold', 
#                       folds = 5, 
#                       nrepeat = n_repeats, 
#                       progressBar = FALSE,
#                       dist = diablo_dist,
#                       auc = TRUE)
#     aucs <- do.call(rbind, diablo_cv$auc)
#     aucs <- data.frame(aucs) %>% 
#       mutate(diablo_design = diablo_design) %>%
#       tibble::rowid_to_column("component") %>%
#       left_join(cors, by = "component")
#     
#     
#     return(aucs)
#   })
#   auc_per_design <- bind_rows(auc_per_design)
#   auc_per_design$signif <- get_signif_mark(auc_per_design$p.value)
#   auc_per_design$component2 <- auc_per_design$component %>% as.character()
#   log_debug("Collected AUC from all designs")
#   
#   comp_colors <- c(c("darkolivegreen2", "darkolivegreen3", "darkolivegreen4", "darkolivegreen", "palegreen4"))[1:diablo_prelim_ncomp]
#   
#   # Plot correlation vs. discrimination tradeoff
#   p1 <- ggplot(auc_per_design, aes(x = diablo_design, y = AUC, group = component2, color = component2)) +
#     geom_line() +
#     geom_point(size = 4) +
#     geom_text(aes(label = signif), color = "black", vjust = -0.3, size = 4) +
#     scale_x_continuous(breaks = designs_to_try) +
#     xlab("") +
#     ylab("DIABLO model AUC") +
#     ggtitle("Model AUC as a function of DIABLO design") +
#     scale_color_manual(values = comp_colors, name = "Number of\ncomponents\nin model") +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5))
#   
#   p2 <- ggplot(auc_per_design, aes(x = diablo_design, y = mean_pairwise_cor, group = component2, color = component2)) +
#     geom_line(position = position_dodge(width=0.06)) +
#     geom_point(size = 4, position = position_dodge(width=0.06)) +
#     geom_errorbar(aes(ymin = mean_pairwise_cor - sd_pairwise_cor, ymax = mean_pairwise_cor + sd_pairwise_cor), width = 0.1, position = position_dodge(width=0.06)) +
#     scale_x_continuous(breaks = designs_to_try) +
#     xlab("DIABLO design constant") +
#     ylab("Mean pairwise correlation") +
#     ggtitle("Mean pairwise correlations between variates") +
#     scale_color_manual(values = comp_colors, name = "Component ID") +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5))
#   
#   p <- plot_grid(p1, p2, ncol = 1, align = "v")
#   log_debug("Done")
#   
#   return(list(p = p, auc_per_design = auc_per_design))
# }
# 
# tune_penalties <- function (diablo_input, diablo_prelim_ncomp, diablo_design, diablo_dist, n_repeats = 20, BPPARAM = SerialParam(), spars_cutoffs = c(4,8,12)) {
#   sparsity_cutoffs <- lapply(diablo_input$X, function(x) spars_cutoffs)
#   
#   diablo_tuning <- tune.block.splsda(
#     X = diablo_input$X, 
#     Y = diablo_input$Y, 
#     ncomp = diablo_prelim_ncomp,
#     test.keepX = sparsity_cutoffs,
#     folds = 10,
#     nrepeat = n_repeats,
#     dist = diablo_dist,
#     design = diablo_design,
#     max.iter = 300,
#     BPPARAM = BPPARAM
#   )
#   
#   # Exploration of selected parameters and their runner-ups
#   tmp_colnames <- paste0(names(diablo_input$X), "_n_vars")
#   tuning_errors_sd <- diablo_tuning$error.rate.sd %>% 
#     data.frame() %>%
#     tibble::rownames_to_column("n_vars_per_omic") %>%
#     tidyr::pivot_longer(
#       cols = colnames(diablo_tuning$error.rate),
#       values_to = "error_rate_sd", names_to = "component"
#       ) 
#   winners <- data.frame(diablo_tuning$choice.keepX)
#   names(winners) <- tmp_colnames
#   winners$final_choice = TRUE
#   winners$component <- paste0('comp', 1:nrow(winners))
#   tuning_errors <- diablo_tuning$error.rate %>% 
#     data.frame() %>%
#     tibble::rownames_to_column("n_vars_per_omic") %>%
#     tidyr::pivot_longer(
#       cols = colnames(diablo_tuning$error.rate),
#       values_to = "error_rate", names_to = "component"
#     ) %>%
#     tidyr::separate(
#       col = "n_vars_per_omic", 
#       sep = "_", 
#       remove = FALSE,
#       into = tmp_colnames,
#       convert = TRUE
#     ) %>%
#     rowwise() %>%
#     mutate(total_vars = sum(!!!syms(tmp_colnames))) %>%
#     left_join(tuning_errors_sd, by = c("n_vars_per_omic", "component")) %>%
#     left_join(winners, by = c("component", tmp_colnames)) %>%
#     mutate(final_choice = tidyr::replace_na(final_choice, FALSE)) %>%
#     arrange(component, error_rate)
#   
#   return(list(tuning_errors = tuning_errors, keepX = diablo_tuning$choice.keepX))
# }
# 
# # Prints tuned parameters for different random seeds. How different are the results?
# # tune_penalties_stability_analysis <- function(diablo_input, diablo_prelim_ncomp, diablo_design, diablo_dist, n_repeats = 3, BPPARAM = SerialParam(), n_seeds = 5) {
# #   log_debug("Analyzing stability of tuned parameters with different random seeds")
# #   for (i in 1:n_seeds) {
# #     log_debug("Tuning penalty parameters with seed {i}")
# #     set.seed(i)
# #     tmp <- tune_penalties(diablo_input, diablo_prelim_ncomp, diablo_design, diablo_dist, n_repeats, BPPARAM)$keepX
# #     log_debug(tmp)
# #   }
# # }
# 
# remove_samples_from_input <- function(diablo_input, held_out_samples) {
#   new_diablo_input <- list()
#   new_diablo_input$X <- lapply(diablo_input$X, function(m) {m[-held_out_samples,]})
#   new_diablo_input$Y <- diablo_input$Y[-held_out_samples]
#   new_diablo_input$sample_ids <- diablo_input$sample_ids[-held_out_samples]
#   return(new_diablo_input)
# }
# 
# get_auc <- function(pred_col, data, label_col = "true_label") {
#   require(pROC)
#   if (n_distinct(data[[label_col]]) == 1) return(NA)
#   f <- as.formula(paste(label_col, "~", pred_col))
#   return(as.numeric(auc(roc(formula = f, data = data, levels = c("healthy", "disease"), direction = "<"))))
# }
# 
# get_cv_auc <- function(diablo_input, rf_cv_results = NULL, folds, list_keepX, final_ncomp, diablo_dist, diablo_design) {
#   require(pROC)
#   require(rsample)
#   
#   n_samples <- length(diablo_input$sample_ids)
#   
#   # Collect AUC stats
#   my_aucs <- data.frame()
#   
#   for (i in 1:nrow(folds)) {
#     fold_id <- paste(folds$splits[[i]]$id, collapse = " ")
#     cat('.')
#     
#     # Partition the data into train and test ("held-out")
#     # held_out_samples_i <- sample(1:n_samples, round(n_samples/10))
#     # train_samples <- setdiff(1:n_samples, held_out_samples_i)
#     
#     train_samples <- folds$splits[[i]]$in_id
#     held_out_samples_i <- setdiff(1:n_samples, train_samples)
#     new_diablo_input <- remove_samples_from_input(diablo_input, held_out_samples_i)
#     held_out_data <- remove_samples_from_input(diablo_input, train_samples)
#     
#     # Train DIABLO model
#     tmp_diablo_out <- block.splsda(
#       new_diablo_input$X, 
#       new_diablo_input$Y, 
#       ncomp = final_ncomp, 
#       keepX = list_keepX, 
#       max.iter = 300,
#       design = diablo_design
#     )
#     
#     # Make prediction on held out samples
#     predictions <- predict(tmp_diablo_out, held_out_data$X, dist = diablo_dist)
#     
#     # Calculate auroc
#     tmp <- data.frame(preds = predictions$AveragedPredict[,1,], true_label = held_out_data$Y)
#     #tmp2 <- bind_rows(tmp2, tmp %>% mutate(fold_id=fold_id))
#     my_aucs_v <- lapply(1:final_ncomp, function(i) get_auc(paste0("preds.dim", i), tmp))
#     names(my_aucs_v) <- paste0("auc_comp", 1:final_ncomp)
#     my_aucs <- bind_rows(
#       my_aucs, 
#       data.frame(my_aucs_v) %>%
#         mutate(fold_id = fold_id) %>%
#         mutate(h_ratio_train = round(unname(table(new_diablo_input$Y)['healthy']) / length(new_diablo_input$Y), 3)) %>%
#         mutate(h_ratio_test = round(unname(table(held_out_data$Y)['healthy']) / length(held_out_data$Y), 3))
#     )
#   }
#   cat('\n')
#   
#   # Add information from ML pipeline results (matched by fold)
#   if (!is.null(rf_cv_results))
#     my_aucs <- my_aucs %>%
#       left_join(
#         rf_cv_results %>% 
#           select(fold_id, out_of_fold_test_auc, mean_out_of_fold_test_auc, h_ratio_train, h_ratio_test) %>% 
#           rename(rf_fold_auc = out_of_fold_test_auc) %>% 
#           rename(mean_rf_auc = mean_out_of_fold_test_auc) %>% 
#           rename(h_ratio_train_SANITY = h_ratio_train, h_ratio_test_SANITY = h_ratio_test),
#         by = "fold_id"
#       )
#   
#   return(my_aucs)
# }
# 
# get_cor_btwn_variates <- function(diablo_out, block1, block2, comp_id) {
#   comp_txts <- paste0('comp', comp_id)
#   if (length(comp_txts) == 1) 
#     comp_txts <- c(comp_txts)
#   unname(sapply(
#     comp_txts, 
#     function(comp_txt) {
#       round(cor(diablo_out$variates[[block1]][,comp_txt], 
#                 diablo_out$variates[[block2]][,comp_txt]), 3) 
#   }))
# }
# 
# get_ttest_variate_vs_label <- function(diablo_out, block1, comp_id) {
#   comp_txts <- paste0('comp', comp_id)
#   if (length(comp_txts) == 1) comp_txts <- c(comp_txts)
#   unname(sapply(comp_txts, function(comp_txt) {
#     variate <- diablo_out$variates[[block1]][,comp_txt]
#     t.test(variate[diablo_out$Y == 'healthy'], variate[diablo_out$Y != 'healthy'])$p.value 
#   }))
# }

# Re-organizes the loadings from a diablo model object into a single compact table
organize_diablo_loadings <- function(diablo_out) {
  require(dplyr)
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
# 
# 
# 
# ###################################
# # Plotting functions
# ###################################
# 
# plot_comp_bubbles <- function(tmp, comp_id) {
#   tmp <- tmp %>%
#     rename_with(function(s) gsub(paste0('comp',comp_id,'_'), '', s))
#   
#   p <- ggplot(tmp, aes(x = '', y = feature_orig2, size = abs(loading), alpha = mean_cv_stability)) +
#     geom_point(color = 'black', fill = "darkslateblue", shape = 21) +
#     ylab(NULL) +
#     #xlab(paste0('DIABLO\ncomponent\n', comp_id)) +
#     xlab(comp_id) +
#     scale_x_discrete(position = "top", expand = expansion(add = .4)) +
#     scale_alpha_continuous(range = c(0.5, 0.8), na.value = 0) +
#     scale_size_continuous(range = c(1.5,4.5)) +
#     # scale_fill_manual(values = c("healthy" = "dodgerblue4", "disease" = "firebrick"),
#     #                   name = 'Increased in...') +
#     theme_minimal() +
#     facet_grid(feat_type ~ ., scales = "free_y", space = "free") +
#     theme(strip.background = element_blank(), 
#           strip.text.x = element_blank(), 
#           strip.text = element_blank()) +
#     theme(legend.position = 'none') +
#     theme(axis.text.y = element_blank()) +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#     theme(panel.grid.major.x = element_line()) +
#     theme(axis.title.x.top = element_text(size = 10, vjust = -2, margin=margin(0,0,0,0)))
#   
#   return(p)
# }
# 
# plot_feature_imp_with_components <- function(
#     d, 
#     diablo_loadings,
#     rf_full_feat_imp,
#     n_comp,
#     # As this graph can be quite long, we limit the number of features 
#     #  we include in the visualization using the following cutoffs:
#     # Top features from each omic that we want to include in visualization
#     rank_cutoff_single_omic = 3, 
#     # Top features from multi-omic models we want to include
#     rank_cutoff_multi_omic = 8, 
#     # This only effects which features from diablo components MUST appear in the list. 
#     #  If other features are already presented and they take part in a component - they will appear as well  
#     stability_cutoff_key_diablo_features = 0.9 
#   ) {
#   
#   # Get feature importance stats
#   tmp <- rf_full_feat_imp %>%
#     filter(dataset == d) %>%
#     select(-dataset)
#   
#   if (nrow(tmp) == 0) {
#     message('Missing information about RF feature importance')
#     return()
#   }
#   
#   # Organize component info
#   tmp2 <- diablo_loadings %>%
#     filter(dataset == d) %>%
#     filter(component <= min(3, n_comp[d])) %>%
#     filter(mean_cv_stability >= 0.1) %>%
#     filter(abs(loading) >= 0.1) %>%
#     mutate(diablo_feature_cluster = ifelse(grepl('_Cluster[0-9]+$', feature),
#                                            gsub("^.*_Cluster([0-9]+)$","\\1",feature),
#                                            NA)) %>%
#     mutate(feature_orig = gsub("_Cluster[0-9]+$","",feature)) # Shortcut...
#   
#   # Patch to mark features that are cluster representatives
#   tmp <- tmp %>%
#     full_join(tmp2 %>% select(feature_orig, diablo_feature_cluster),
#               by = "feature_orig")
#   
#   # Add component info (stable components and features only)
#   for (i in 1:n_comp[d]) {
#     tmp_comp <- tmp2 %>%
#       filter(component == i) %>%
#       select(-dataset, -feature, -component, -feature_set, -diablo_feature_cluster) %>%
#       rename_with(.fn = function(s) paste0('comp',i,'_',s), 
#                   .cols = all_of(c("loading", "mean_cv_stability")))
#     
#     if (nrow(tmp_comp) == 0) next;
#     
#     tmp <- tmp %>%
#       full_join(tmp_comp, by = c("feature_orig"))
#   }
#   
#   # Add column to indicate whether a feature participates in DIABLO components
#   tmp <- tmp %>%
#     mutate(key_in_diablo_comps = feature_orig %in% 
#              (tmp2 %>% filter(mean_cv_stability > stability_cutoff_key_diablo_features) %>% pull(feature_orig)))
#   
#   # Keep only features either ranked highly in RF models or taking part in diablo components
#   tmp <- tmp %>%
#     # Top X from each pipeline
#     filter(single_omic_rank <= rank_cutoff_single_omic | 
#              multi_omic_rank <= rank_cutoff_multi_omic |
#              key_in_diablo_comps) %>%
#     # Shorten names of features that are too long
#     mutate(feature_orig2 = ifelse(nchar(feature_orig) > 40,
#                                   paste0(substr(feature_orig, 1, 40),"..."),
#                                   feature_orig)) %>%
#     # Mark cluster representatives (in either pipeline) with a star
#     mutate(
#       feature_orig2 = 
#         paste0(
#           feature_orig2,
#           ifelse((!is.na(feature_pretty_single_omic) & 
#                     grepl('_Cluster[0-9]+$', feature_pretty_single_omic)) |
#                    (!is.na(feature_pretty_multi_omic) & 
#                       grepl('_Cluster[0-9]+$', feature_pretty_multi_omic)) |
#                    (!is.na(diablo_feature_cluster)),
#                  '*', ''))) %>%
#     # Once marked we erase unneeded columns
#     select(-diablo_feature_cluster) %>%
#     mutate(feat_type = substr(feature_orig2,1,1)) %>%
#     mutate(feat_type = factor(feat_type, levels = c('T','G','P','M'))) 
#   
#   # Re-run u-tests per feature. TODO: fdr should be applied before subsetting to only these top features
#   tmp$increased_in <- NA
#   tmp$fdr_utest <- NA
#   
#   for (i in 1:nrow(tmp)) {
#     # Get feature name and type
#     f <- tmp$feature_orig[i]
#     f_type <- as.character(tmp$feat_type)[i]
#     
#     # Extract actual data
#     df <- diablo_out[[d]][['X']][[f_type]]
#     
#     ##########################################################################
#     # PATCH! 
#     # TODO: f won't be found in df (below) in cases where the diablo clusters 
#     #  were a bit different from the clusters in the other pipelines...
#     #  Need to read diablo clusters file and map f to cluster representative.
#     if (d == "crc_feng_2015" & f == "P__GLUTORN.PWY") f <- "P__ARGSYNBSUB.PWY"
#     if (d == "crc_feng_2015" & f == "P__PWY.4242") f <- "P__COA.PWY"
#     if (d == "cdi_schubert_16s_2014" & f == "P__MET.SAM.PWY") f <- "P__HOMOSER.METSYN.PWY"
#     if (d == "sth_rubel_2020" & f == "M__C04932") f <- "M__C04919"
#     print(f)
#     ##########################################################################
#     
#     # Remove cluster identifiers from column names (just for the sake of easily extracting the feature)
#     colnames(df) <- gsub('_Cluster[0-9]+$','',colnames(df))
#     healthy_ind <- diablo_out[[d]]$Y == "healthy"
#     h_vec <- df[,f][healthy_ind]
#     d_vec <- df[,f][!healthy_ind]
#     
#     # Run the u-tests (TODO: replace)
#     utest <- wilcox.test(h_vec, d_vec, digits.rank = 6) 
#     tmp$fdr_utest[i] <- utest$p.value
#     
#     # If test is significant, mark which group has increased values
#     if (utest$p.value <= 0.1) {
#       tmp$increased_in[i] <- ifelse(mean(h_vec) > mean(d_vec),
#                                     "Healthy", "Disease")
#     } else {
#       tmp$increased_in[i] <- "Non-significant"
#     }
#   }
#   
#   # Correct p values
#   tmp <- tmp %>%
#     mutate(fdr_utest = p.adjust(fdr_utest, method = "fdr")) %>%
#     mutate(fdr_utest_for_plot = -log10(fdr_utest)) 
#   
#   # Plot 1: increased in healthy or disease
#   p1 <- ggplot(tmp, aes(x = '', y = feature_orig2, fill = increased_in)) +
#     geom_point(color = 'black', shape = 21, size = 4) +
#     ylab(NULL) +
#     xlab(NULL) +
#     scale_x_discrete(position = "top", expand = expansion(add = .4)) +
#     scale_fill_manual(values = c("Healthy" = "dodgerblue4", 
#                                  "Disease" = "firebrick", 
#                                  "Non-significant" = "lightgrey"),
#                       name = 'Increased in...') +
#     theme_minimal() +
#     facet_grid(feat_type ~ ., scales = "free_y", space = "free") +
#     theme(strip.background = element_blank(), 
#           strip.text.x = element_blank(), 
#           strip.text = element_blank()) +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#   
#   # Plot 2
#   p2 <- ggplot(tmp, aes(x = '', y = feature_orig2, fill = fdr_utest_for_plot)) +
#     geom_tile(color = 'black') +
#     ylab(NULL) +
#     xlab('Univariate\ntest') +
#     scale_x_discrete(position = "top", expand = c(0, 0)) +
#     scale_fill_gradient(high = 'deepskyblue4', 
#                         low = lighten('aliceblue'),
#                         name = 'U-test FDR\n(-log10)',
#                         guide = guide_colorbar(frame.colour = "black", 
#                                                ticks.colour = "black")) +
#     theme_minimal() +
#     facet_grid(feat_type~., scales = "free_y", space = "free") +
#     theme(strip.background = element_blank(), 
#           strip.text.x = element_blank(), 
#           strip.text = element_blank()) +
#     theme(axis.text.y = element_blank()) +
#     theme(axis.title.x.top = element_text(size = 10, vjust = -2, margin=margin(0,0,0,0)))
#   
#   # Plot 3
#   p3 <- ggplot(tmp, aes(x = '', y = feature_orig2, 
#                         fill = single_omic_mean_importance)) +
#     geom_tile(color = 'black') +
#     ylab(NULL) +
#     xlab('Single-omic\nfeature\nimportance') +
#     scale_x_discrete(position = "top", expand = c(0, 0)) +
#     scale_fill_gradient(high = 'darkolivegreen', 
#                         low = lighten('darkolivegreen1', amount = 0.8),
#                         na.value = "grey85",
#                         name = 'Feature importance\nscore -\nSingle-omic model',
#                         guide = guide_colorbar(frame.colour = "black", 
#                                                ticks.colour = "black")) +
#     theme_minimal() +
#     facet_grid(feat_type~., scales = "free_y", space = "free") +
#     theme(strip.background = element_blank(), 
#           strip.text.x = element_blank(), 
#           strip.text = element_blank()) +
#     theme(axis.text.y = element_blank()) +
#     theme(axis.title.x.top = element_text(size = 10, vjust = -2, margin=margin(0,0,0,0)))
#   
#   # Plot 4
#   p4 <- ggplot(tmp, aes(x = '', y = feature_orig2, fill = multi_omic_mean_importance)) +
#     geom_tile(color = 'black') +
#     ylab(NULL) +
#     xlab('Multi-omic\nfeature\nimportance') +
#     scale_x_discrete(position = "top", expand = c(0, 0)) +
#     scale_fill_gradient(high = 'darkorange4', low = 'blanchedalmond',
#                         na.value = "grey85",
#                         name = 'Feature importance\nscore -\nMulti-omic model',
#                         guide = guide_colorbar(frame.colour = "black", 
#                                                ticks.colour = "black")) +
#     theme_minimal() +
#     facet_grid(feat_type~., scales = "free_y", space = "free") +
#     theme(strip.background = element_blank(), 
#           strip.text.x = element_blank(), 
#           strip.text = element_blank()) +
#     theme(axis.text.y = element_blank()) +
#     theme(axis.title.x.top = 
#             element_text(size = 10, vjust = -2, margin=margin(0,0,0,0)))
#   
#   p5 <- plot_comp_bubbles(tmp, 1)
#   if (n_comp[d]>1) p6 <- plot_comp_bubbles(tmp, 2) else p6 <- ggplot() + theme_minimal()
#   if (n_comp[d]>2) p7 <- plot_comp_bubbles(tmp, 3) else p7 <- ggplot() + theme_minimal()
#   
#   # Extract the legends
#   p_legends <- plot_grid(
#     get_legend(p1 + theme(legend.title = element_text(size=10), legend.box.margin = margin(10, 0, 0, 15))),
#     get_legend(p2 + theme(legend.title = element_text(size=10), legend.box.margin = margin(0, 0, 0, 15))),
#     get_legend(p3 + theme(legend.title = element_text(size=10), legend.box.margin = margin(0, 0, 0, 15))),
#     get_legend(p4 + theme(legend.title = element_text(size=10), legend.box.margin = margin(0, 0, 0, 15))),
#     ncol = 1,
#     align = 'v',
#     rel_heights = c(1,2,2,2)
#   )
#   
#   p <- plot_grid(p1 + theme(legend.position = "none"), 
#                  p2 + theme(legend.position = "none"), 
#                  p3 + theme(legend.position = "none"), 
#                  p4 + theme(legend.position = "none"), 
#                  p5,
#                  p6,
#                  p7,
#                  p_legends,
#                  rel_widths = c(6,1.5,1.5,1.5,0.7,0.7,0.7,3), 
#                  nrow = 1,
#                  align = 'h', 
#                  axis = 'tb')
#   
#   print(p)
#   return(tmp)
# }

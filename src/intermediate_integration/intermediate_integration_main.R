###################################
# Read command line arguments
###################################

library(optparse)
library(BiocParallel)
library(logger)

option_list = list(
  make_option(
    c("-d", "--dataset"), 
    type = "character", 
    default = "cd_franzosa_2019", 
    help = "dataset name", 
    metavar = "character"
  ),
  make_option(
    c("-w", "--work_dir"), 
    type = "character", 
    default = getwd(), 
    help = "working directory (should be set to the ml_pipeline directory)", 
    metavar = "character"
  )
); 

opt <- parse_args(OptionParser(option_list=option_list))
d <- opt$dataset
log_debug("Starting intermediate integration analysis of dataset {d}")
setwd(opt$work_dir)

###################################
# Preparations
###################################

library(config)
library(mixOmics)
library(dplyr)
library(readr)
library(rsample)
library(pROC)
library(igraph)
library(ranger)

source('clustering.R')
source("preprocessing.R")
source("../diablo/diablo_utils.R")
log_threshold('DEBUG')
options(scipen = 999)

###################################
# Define output files 
###################################

res_dir <- '../../data/intermediate_integration_results'
dir.create(res_dir, showWarnings = FALSE)
out_clusters <- paste0(res_dir, '/', d, "_clusters.csv") # Ignore
out_rdata <- paste0(res_dir, '/', d, '_R_objects.RData')

###################################
# Load and organize data 
###################################

# Load processed data
log_debug("Loading dataset: {d}")
include_metabs <- TRUE

# File names
metadata_path <- sprintf(config::get("paths_templates")$metadata, d)
taxonomy_path <- sprintf(config::get("paths_templates")$taxonomy, d)
pathways_path <- sprintf(config::get("paths_templates")$pathways, d)
if (include_metabs) metabolites_path <- sprintf(config::get("paths_templates")$metabolites, d)
metagenome_path <- sprintf(config::get("paths_templates")$metagenome, d)

# Load datasets and remove rare/constant features
metadata_df <- prep_load_metadata(metadata_path) # Only two columns retained: 'sample_id__' and 'DiseaseState'
taxonomy_df <- prep_sanitize_dataset(prep_load_taxonomy(taxonomy_path),d) 
pathways_df <- prep_sanitize_dataset(prep_load_pathways(pathways_path),d) 
if (include_metabs) metabolites_df <- prep_sanitize_dataset(prep_load_metabolites(metabolites_path),d) 
genes_df <- prep_sanitize_dataset(prep_load_metagenome(metagenome_path),d)

# Fix names
colnames(taxonomy_df) <- make.names(colnames(taxonomy_df))
colnames(pathways_df) <- make.names(colnames(pathways_df))
if (include_metabs) colnames(metabolites_df) <- make.names(colnames(metabolites_df))
colnames(genes_df) <- make.names(colnames(genes_df))

# CLR transpose relevant tables
ps_count <- 0.0000001
taxonomy_df_clr <- prep_clr_transpose(taxonomy_df, ps_count)
genes_df_clr <- prep_clr_transpose(genes_df, ps_count)
pathways_df_clr <- prep_clr_transpose(pathways_df, ps_count)

# Log transform metabolites
# (First - replace zeros/na's with half the minimal non-zero value)
if (include_metabs) 
  metabolites_df_log <- metabolites_df %>%
    mutate_if(is.numeric, function(x) ifelse(x==0, NA, x)) %>%
    mutate_if(is.numeric, function(x) ifelse(is.na(x), min(x, na.rm = T)/2, x)) %>%
    mutate_if(is.numeric, log)

# Join 
if (include_metabs) {
  tmp <- metadata_df %>%
    inner_join(taxonomy_df_clr, by = "sample_id__") %>%
    inner_join(pathways_df_clr, by = "sample_id__") %>%
    inner_join(metabolites_df_log, by = "sample_id__") %>%
    inner_join(genes_df_clr, by = "sample_id__") 
} else {
  tmp <- metadata_df %>%
    inner_join(taxonomy_df_clr, by = "sample_id__") %>%
    inner_join(pathways_df_clr, by = "sample_id__") %>%
    inner_join(genes_df_clr, by = "sample_id__")
}

# Remove highly-redundant features and replace with cluster representatives,
#  as previously calculated
tmp2 <- cluster_cluster_features(
  dataset = tmp %>% select(-sample_id__), 
  feature_set_type = ifelse(include_metabs, "T+G+P+M", "T+G+P"),
  cluster_type = "clustering99",
  clusters_output = out_clusters
)
proc_data <- cbind(tmp %>% select(sample_id__), tmp2)
log_debug("Data loaded")

###################################
# Intermediate integration
###################################

# 1. Organize input to DIABLO
# ------------------------------------------------
diablo_input <- organize_data_for_diablo2(d, proc_data)
log_debug("Completed data preparations for DIABLO")
# summary(diablo_input$Y)

# 2. Wrap DIABLO in a repeated sub-sampling procedure 
#  to identify robust components
# ------------------------------------------------

log_debug('Running sensitivity analysis over component identification parameters')

# Function for wrapped DIABLO (using specific pipeline parameters)
run_wrapped_diablo <- function(
    diablo_input,
    diablo_design = 2/3,
    n_repeats = 3,
    n_folds = 10,
    diablo_keepX = 10, 
    diablo_ncomp = 5,
    diablo_dist = 'centroids.dist') {
  
  # Set number of variables to take from each omic
  list_keepX <- list()
  for(x in names(diablo_input$X)) list_keepX[[x]] <- rep(diablo_keepX, diablo_ncomp)
  
  # Get number of samples
  n_samples <- length(diablo_input$sample_ids)
  
  # Generate data subsets (=repeated cross-validation) 
  set.seed(27) 
  folds <- vfold_cv(
    data.frame(sample_num = 1:n_samples), 
    v = n_folds, 
    repeats = n_repeats
  )
  
  # Place holders for AUC stats and component loadings
  cv_loadings <- data.frame()
  
  for (i in 1:nrow(folds)) {
    fold_id <- paste(folds$splits[[i]]$id, collapse = "_")
    cat('.')
    
    # Partition the data into train and test ("held-out")
    train_samples <- folds$splits[[i]]$in_id
    held_out_samples_i <- setdiff(1:n_samples, train_samples)
    new_diablo_input <- remove_samples_from_input(diablo_input, held_out_samples_i)
    held_out_data <- remove_samples_from_input(diablo_input, train_samples)
    
    # Train DIABLO model
    tmp_diablo_out <- block.splsda(
      new_diablo_input$X, 
      new_diablo_input$Y, 
      ncomp = diablo_ncomp, 
      keepX = list_keepX, 
      max.iter = 300,
      design = diablo_design
    )
    
    # Make prediction on held out samples using DIABLO's approach 
    # predictions <- predict(tmp_diablo_out, held_out_data$X, dist = diablo_dist)
    
    # Calculate auroc on held out data
    # tmp <- data.frame(
    #   preds = predictions$AveragedPredict[,1,], 
    #   true_label = held_out_data$Y
    #   )
    # cv_aucs_v <- lapply(
    #   1:diablo_ncomp, 
    #   function(i) get_auc(paste0("preds.dim", i), tmp)
    #   )
    # names(cv_aucs_v) <- paste0("auc_comp", 1:diablo_ncomp)
    # 
    # cv_aucs <- bind_rows(
    #   cv_aucs, 
    #   data.frame(cv_aucs_v) %>%
    #     mutate(fold_id = fold_id) 
    # )
    
    # Save all non-zero loadings
    tmp_loadings <- organize_diablo_loadings(tmp_diablo_out) %>%
      # Add full component identifiers (i.e. component number + specific fold)
      mutate(fold_id = fold_id) %>%
      mutate(comp_fold_id = paste0(fold_id, '_Comp', component)) %>%
      dplyr::select(-component)
    cv_loadings <- bind_rows(cv_loadings, tmp_loadings)
  }
  cat('\n')
  
  # For each pair of features count the number of times they appeared together
  # (Note - each pair currently appears twice: [A,B], [B,A])
  feat_pairs <- 
    # Get pairs of features that were found in the same component
    cv_loadings %>%
      full_join(cv_loadings, by = c("fold_id","comp_fold_id")) %>%
      filter(feature.x != feature.y) %>%
      # In cases where a pair of features come out together in more than one 
      #  component in the same run, we only count it once
      dplyr::select(feature.x, feature_set.x, feature.y, feature_set.y, fold_id) %>% 
      distinct() %>%
      group_by(feature.x, feature.y) %>%
      summarise(N = n(), .groups = "drop") %>%
      # Compute a distance between each pair, where 1 = never appeared together, 
      #  i.e. maximal distance, 0 = always appear together, minimal distance
      mutate(dist = ((n_repeats * n_folds) - N) / (n_repeats * n_folds)) 
  
  return(list(feat_pairs = feat_pairs))
}

# Run sensitivity analysis on:
# - keepX
# - design
# - n_repeats/folds
# - clustering threshold
feat_pairs <- data.frame()

# For convenience - create a string that represents the entire parameter choices
run_id_pattern <- 'keep_%i//des_%g//nrep_%i//nfol_%i//ncom_%i//edge_%g'

for (diablo_keepX in c(5, 10, 15)) {
  for (diablo_design in c(0.666, 1)) {
    for (n_repeats in c(10)) {
      for (n_folds in c(5, 10)) {
        for (diablo_ncomp in c(5)) {
          print(paste(
            diablo_keepX,
            diablo_design,
            n_repeats,
            n_folds,
            diablo_ncomp,
            sep = ' // '
          ))
          
          # Run wrapped diablo
          tmp <- run_wrapped_diablo(
            diablo_input,
            diablo_design = diablo_design,
            n_repeats = n_repeats,
            n_folds = n_folds,
            diablo_keepX = diablo_keepX, 
            diablo_ncomp = 10
          )
          
          feat_pairs <- bind_rows(
            feat_pairs,
            tmp$feat_pairs %>%
              mutate(param_keepX = diablo_keepX) %>%
              mutate(param_diablo_design = diablo_design) %>%
              mutate(param_n_repeats = n_repeats) %>%
              mutate(param_n_folds = n_folds) %>%
              mutate(param_ncomp = diablo_ncomp) 
          )
        } # Done iterating over diablo_ncomp
      } # Done iterating over n_folds
    } # Done iterating over n_folds
  } # Done iterating over n_folds
} # Done iterating over diablo_keepX

# 3. Extract final components for each setting 
#  (for sensitivity analysis)
# ------------------------------------------------

# We add one last parameter for our sensitivity analysis:
# the 'edge_threshold' defines for a pair of 2 features, what portion of 
# sub-sampling repeats do we require these features to be included in the same 
# module in order for them to be placed together in our final modules.
log_debug('Selecting final compoments')
tmp_feat_pairs <- bind_rows(
  feat_pairs %>% mutate(param_edge_threshold = 0.7),
  feat_pairs %>% mutate(param_edge_threshold = 0.8)
  ) %>%
  # Add a "run identifier"
  mutate(run_id = sprintf(
    run_id_pattern, 
    param_keepX, 
    param_diablo_design, 
    param_n_repeats, 
    param_n_folds, 
    param_ncomp, 
    param_edge_threshold
  ))
# table(tmp_feat_pairs$run_id)

# Thin table to record all runs for sensitivity analysis
sens_analysis_runs <- tmp_feat_pairs %>%
  select(run_id, starts_with('param_')) %>%
  distinct() %>%
  tidyr::pivot_longer(
    cols = starts_with('param_'), 
    names_to = 'param', 
    names_prefix = 'param_', 
    values_to = 'value'
    )

# For each setting, infer modules and organize them in a list
sens_analysis_components <- data.frame()
for (curr_set in unique(sens_analysis_runs$run_id)) {
  cat('.')
  # Get number of overall repeats
  n_rep_folds <- sens_analysis_runs %>%
    filter(run_id == curr_set) %>%
    filter(param %in% c('n_folds','n_repeats')) %>%
    pull(value)
  n_rep_folds <- n_rep_folds[1] * n_rep_folds[2]
  
  # Get edge threshold
  edge_threshold <- sens_analysis_runs %>%
    filter(run_id == curr_set) %>%
    filter(param == 'edge_threshold') %>%
    pull(value)
    
  # Create list of edges by taking all feature pairs that appeared together 
  #  over <edge_threshold> of the times
  edges <- tmp_feat_pairs %>% 
    filter(run_id == curr_set) %>%
    filter(N >= n_rep_folds * edge_threshold)
  
  if(nrow(edges) == 0) {
    cat('\n')
    log_debug('No components identified with setting {curr_set}')
    next
  }
  
  vertices <- edges %>% 
    select(feature.x) %>% 
    distinct() 
  
  net <- graph_from_data_frame(
    d = edges, 
    vertices = vertices, 
    directed = F
  ) 
  
  tmp <- components(net, mode = "strong")
  final_components <- data.frame(
    feature = names(tmp$membership),
    component = paste0('comp',unname(tmp$membership)),
    run_id = curr_set
  ) %>% tibble::remove_rownames()
  
  sens_analysis_components <- bind_rows(
    sens_analysis_components,
    final_components
  )
}
cat('\n')
rm(tmp, final_components, net, vertices, edges, curr_set)

# Discard settings for which no modules were identified
sens_analysis_runs <- sens_analysis_runs %>%
  filter(run_id %in% unique(sens_analysis_components$run_id))

# 4. Compute module summaries (1st PC), per setting
# --------------------------------------------------

log_debug('Computing new variables based on components, and shuffled versions for comparison')
latent_vars <- list()
latent_var_t_tests <- list()
comp_variance_expl <- list()
comp_inter_omic_cors <- list()
shuffled_components <- list()
set.seed(27)

# For each setting, get modules first PCs
for (curr_set in unique(sens_analysis_runs$run_id)) {
  
  # Initialize
  latent_vars[[curr_set]] <- list()
  latent_var_t_tests[[curr_set]] <- list()
  comp_variance_expl[[curr_set]] <- list()
  comp_inter_omic_cors[[curr_set]] <- list()
  shuffled_components[[curr_set]] <- list()
  
  # Organize current settings 
  tmp <- strsplit(strsplit(curr_set, "//")[[1]], '_')
  selected_settings <- lapply(tmp, function(x) as.numeric(x[2]))
  names(selected_settings) <- sapply(tmp, function(x) x[1])
  
  # Get actual component list
  selected_components <- sens_analysis_components %>%
    filter(run_id == curr_set) %>%
    select(-run_id)
  log_debug('With settings {curr_set}, identified {n_distinct(selected_components$component)} components, with {selected_components %>% group_by(component) %>% summarise(N = n()) %>% pull(N) %>% mean() %>% round(2)} features on average')
  
  # 0 = true components, 1-99 = shuffled versions
  # log_debug('Fetching shuffled components for null distributions')
  for (i in 0:99) {
    shuf <- ifelse(i == 0, FALSE, TRUE)
    run <- ifelse(shuf, paste0('shuf', i), 'true')
    
    latent_vars[[curr_set]][[run]] <- data.frame(label = diablo_input$Y) # Place holder
    comp_variance_expl[[curr_set]][[run]] <- data.frame()
    comp_inter_omic_cors[[curr_set]][[run]] <- data.frame()
    shuffled_components[[curr_set]][[run]] <- list()
    
    # Iterate over components 
    for (comp_id in 1:n_distinct(selected_components$component)) {
      # message("Component ", comp_id)
      
      # Get component
      tmp <- selected_components %>%
        filter(component == paste0('comp', comp_id))
      
      # Extract relevant data to use as PCA input
      # (take only features in this component from each omic)
      xlist <- lapply(diablo_input$X, function(x) {x[,colnames(x) %in% tmp$feature,drop=FALSE]})
      
      # Remove omics without features in this component
      rem <- sapply(xlist, ncol) == 0
      xlist <- xlist[!rem]
      
      # Shuffle if needed 
      if (shuf) {
        # For each omic, sample same amount of features as in true components
        xlist <- lapply(
          names(diablo_input$X),
          function(fs) {
            num_feats_to_sample <- sum(grepl(paste0("^",fs,"__"), tmp$feature))
            x <- diablo_input$X[[fs]]
            x[,sample(colnames(x), num_feats_to_sample),drop=FALSE]
          }
        )
        names(xlist) <- names(diablo_input$X)
        shuffled_components[[curr_set]][[run]][[comp_id]] <- unname(unlist(lapply(xlist, colnames)))
      }
      
      # Run PCA on features in component
      comp_data <- bind_cols(lapply(xlist, as.data.frame))
      pca <- prcomp(comp_data, center = TRUE, scale. = TRUE) 
      var_explained <- pca$sdev^2/sum(pca$sdev^2) * 100 # Percent of variance explained
      
      # Record variance explained (we expect high values for 1st components)
      comp_variance_expl[[curr_set]][[run]] <- bind_rows(
        comp_variance_expl[[curr_set]][[run]],
        data.frame(comp = comp_id,
                   var_explained_pc1 = var_explained[1],
                   var_explained_pc2 = var_explained[2]) %>% 
          tibble::remove_rownames()
      ) 
      
      # Record mean correlations between features from different omics
      tmp_inter_omic_cors <- get_all_inter_omic_corrs(comp_data)
      comp_inter_omic_cors[[curr_set]][[run]] <- bind_rows(
        comp_inter_omic_cors[[curr_set]][[run]],
        data.frame(comp = comp_id,
                   n_pairs = nrow(tmp_inter_omic_cors),
                   avg_pears_corr = mean(tmp_inter_omic_cors$pearson_corr),
                   avg_spear_corr = mean(tmp_inter_omic_cors$spearman_corr))
      ) 
      
      # Extract and save new latent feature
      new_latent <- data.frame(pca$x[,'PC1'])
      names(new_latent) <- paste0("comp", comp_id)
      latent_vars[[curr_set]][[run]] <- bind_cols(
        latent_vars[[curr_set]][[run]], 
        new_latent
      )
    }
    
    # Get t-test p values for latent variables
    latent_var_t_tests[[curr_set]][[run]] <- sapply(
      colnames(latent_vars[[curr_set]][[run]])[-1],
      function(col)
        t.test(as.formula(paste(col," ~ label")), 
               latent_vars[[curr_set]][[run]])$p.value
    ) %>%
      as.data.frame() %>%
      rename(t_test_p_value = 1) %>%
      tibble::rownames_to_column("variate")
    
  }
  
  # Re-organize
  latent_var_t_tests[[curr_set]] <- bind_rows(latent_var_t_tests[[curr_set]], .id = 'run')
  comp_variance_expl[[curr_set]] <- bind_rows(comp_variance_expl[[curr_set]], .id = 'run')
  comp_inter_omic_cors[[curr_set]] <- bind_rows(comp_inter_omic_cors[[curr_set]], .id = 'run')
}


# 5. Run repeated cross validation to assess AUC
# ------------------------------------------------

# First, get the number of repeats/folds from ML results
log_debug('Starting model evaluation using cross-validation')
n_ml_repeats <- config::get('outer_n_repeats')
n_ml_folds <- config::get('outer_n_folds')
n_samples <- length(diablo_input$sample_ids)

# Now generate folds consistent with ML pipeline
set.seed(1111)
folds <- vfold_cv(
  data.frame(sample_num = 1:n_samples),
  v = n_ml_folds,
  repeats = n_ml_repeats
)

# Collect AUC stats here
cv_aucs <- data.frame()

# Iterate over settings
for (curr_set in unique(sens_analysis_runs$run_id)) {
  log_debug('Training and evaluating models using settings {curr_set}')
  
  # Get actual component list
  selected_components <- sens_analysis_components %>%
    filter(run_id == curr_set) %>%
    select(-run_id)
  
  # Iterate over true components and shuffled ones
  for (run in names(latent_vars[[curr_set]])) {
    cat('.')
    
    # Cross validation loop
    for (i in 1:nrow(folds)) {
      fold_id <- paste(folds$splits[[i]]$id, collapse = " ")
      
      # Get train samples
      train_samples <- folds$splits[[i]]$in_id
      # train_data <- latent_vars[[curr_set]][[run]][train_samples,] %>%
      #   # Important for GLM
      #   mutate(label = factor(label, levels = c("healthy", "disease")))
      # test_data <- latent_vars[[curr_set]][[run]][-train_samples,] %>%
      #   mutate(label = factor(label, levels = c("healthy", "disease")))
      
      train_data <- data.frame(label = diablo_input$Y[train_samples])
      test_data <- data.frame(label = diablo_input$Y[-train_samples])
      
      # Iterate over components
      for (comp_id in 1:n_distinct(selected_components$component)) {
        # Get only features in component
        tmp <- selected_components %>%
          filter(component == paste0('comp', comp_id)) %>%
          pull(feature)
        
        # If we're running on shuffled components - take shuffled component features (override real features)
        if (run != 'true') {
          tmp <- shuffled_components[[curr_set]][[run]][[comp_id]]
        }
        comp_data <- lapply(diablo_input$X, function(x) {x[,colnames(x) %in% tmp, drop=FALSE]})
        
        # Remove omics without features in this component
        rem <- sapply(comp_data, ncol) == 0
        comp_data <- bind_cols(lapply(comp_data[!rem], as.data.frame))
        
        # Train PCA
        pca <- prcomp(comp_data[train_samples,], center = TRUE, scale. = TRUE) 
        
        # Run same PCA loadings on test data (test data gets scaled according to train data)
        predicted_pcs <- predict(pca, newdata = comp_data[-train_samples,])
        
        # Organize train + test tables
        tmp_train_data <- pca$x[,'PC1'] %>% data.frame(); names(tmp_train_data)[1] <- paste0('comp',comp_id)
        train_data <- bind_cols(train_data, tmp_train_data)
        tmp_test_data <- predicted_pcs[,'PC1'] %>% data.frame(); names(tmp_test_data)[1] <- paste0('comp',comp_id)
        test_data <- bind_cols(test_data, tmp_test_data)
      }
      
      # Important for GLM
      train_data <- train_data %>%
        mutate(label = factor(label, levels = c("healthy", "disease")))
      test_data <- test_data %>%
        mutate(label = factor(label, levels = c("healthy", "disease")))
      
      # Train logistic regression / RF
      logistic_model <- glm(label ~ ., data = train_data, family = "binomial")
      rf <- ranger(label ~ ., data = train_data, probability = TRUE)
      
      # Make prediction on held out samples
      preds_glm <- predict(logistic_model, test_data, type = "response")
      preds_rf <- predict(rf, test_data)$predictions
      
      # Calculate and save auroc
      tmp <- data.frame(preds_glm = preds_glm, preds_rf = preds_rf[,'disease'], true_label = test_data$label)
      cv_aucs <- bind_rows(
        cv_aucs,
        data.frame(
          setting = curr_set,
          run = run,
          fold_id = fold_id,
          test_auc_glm = get_auc('preds_glm', tmp),
          test_auc_rf = get_auc('preds_rf', tmp)
        )
      )
    }
  }
  cat('\n')
}
rm(tmp, preds_glm, preds_rf, train_samples, train_data, test_data, logistic_model, rf, fold_id, n_samples, n_ml_folds, n_ml_repeats)

# Summarise
summary_aucs <- cv_aucs %>%
  group_by(setting, run) %>%
  summarise(mean_comp_rf_auc = mean(test_auc_rf, na.rm = TRUE),
            sd_comp_rf_auc = sd(test_auc_rf, na.rm = TRUE),
            mean_comp_glm_auc = mean(test_auc_glm, na.rm = TRUE),
            sd_comp_glm_auc = sd(test_auc_glm, na.rm = TRUE),
            .groups = "drop")

# 9. Organize and save results in an RData file
# ------------------------------------------------

save(diablo_input, 
     feat_pairs,
     sens_analysis_components,
     sens_analysis_runs,
     latent_vars,
     latent_var_t_tests,
     comp_variance_expl,
     comp_inter_omic_cors,
     cv_aucs,
     summary_aucs,
     file = out_rdata)
log_debug('Saved diablo objects: {out_rdata}')


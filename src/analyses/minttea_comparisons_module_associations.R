# udocker run --volume=/specific/elhanan/PROJECTS/MULTI_VIEW_EM/repo/multi_view_microbiome_analysis efrat_ubun_r Rscript /specific/elhanan/PROJECTS/MULTI_VIEW_EM/repo/multi_view_microbiome_analysis/src/diablo/minttea_comparisons_module_associations.R
setwd('/specific/elhanan/PROJECTS/MULTI_VIEW_EM/repo/multi_view_integration_analysis')
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) stop("At least one argument must be supplied (dataset name).", call.=FALSE)
d <- args[1] # d <- 'cd_franzosa_2019'

# Preperations
library(mixOmics)
library(dplyr)
library(logger)
library(PMA)
library(readr)
library(igraph)
library(pROC)
library(ppcor)

log_threshold('DEBUG')
set.seed(2727)
options(scipen = 999)

source("src/analyses/scca_and_diablo_utils.R")
source("src/intermediate_integration/MintTea.R")

# Load some dataset (after processing)
res_dir <- 'data/intermediate_integration_results'
f_path <- file.path(res_dir, paste0(d,'_results.RData')) 
load(f_path, envir = tmp_env <- new.env())
data_X <- tmp_env$minttea_results$diablo_input$X
data_Y <- factor(tmp_env$minttea_results$diablo_input$Y, levels = c('healthy','disease'))
data_sample_ids <- tmp_env$minttea_results$diablo_input$sample_ids
n_samples_tot <- length(data_sample_ids)

# As components get less and less stable, we take only the first 3
n_comps <- 3

# Number of repeats
n_repeats <- 100

# Threads
# n_threads = 5

# Additional parameters
perc_samples_to_exclude <- 30
diablo_design <- .5 # To make the comparison against MultiCCA (that ignores study group) more fair
minttea_edge_threshold <- .7
n_minttea_repeats <- 10

n_excluded <- floor(n_samples_tot * perc_samples_to_exclude / 100)
n_kept <- n_samples_tot - n_excluded

####################################################################
# DIABLO
####################################################################

# Results place holder
component_associations <- data.frame()

# Iterate over a few different diablo keepX params
for (diablo_keepX in c(7, 10)) { # diablo_keepX <- 10 
  log_debug('diablo_keepX = ', diablo_keepX)
  
  # Format to diablo's requirement
  list_keepX <- list()
  for(x in names(data_X)) list_keepX[[x]] <- rep(diablo_keepX, n_comps)
  
  # Run 100 times 
  for (i in 1:n_repeats) {
    
    # Randomly sample samples (rows)
    samples_subset <- sample(x = n_samples_tot, size = n_kept)
    samples_oof <- setdiff(1:n_samples_tot, samples_subset)
    data_X_train <- lapply(data_X, function(m, samps) return(m[samps,]), samps = samples_subset)
    data_Y_train <- data_Y[samples_subset]
    
    data_X_test <- lapply(data_X, function(m, samps) return(m[samps,]), samps = samples_oof)
    data_Y_test <- data_Y[samples_oof]
    
    # Run DIABLO for data subset 1
    tmp_diablo_out <- block.splsda(
      data_X_train, 
      data_Y_train, 
      ncomp = n_comps, 
      keepX = list_keepX, 
      max.iter = 300,
      design = diablo_design
    )
    
    # Predict components for oof data
    predictions <- predict(object = tmp_diablo_out, newdata = data_X_test) 
    predictions <- predictions$variates
    
    # Get loadings
    loadings <- organize_diablo_loadings(tmp_diablo_out) 
    
    # Compute correlations between variates for oof data
    comp_data <- bind_cols(lapply(data_X_test, as.data.frame))
    
    # Record mean correlations between features from different omics
    mean_inter_omic_cors <- c()
    for (j in 1:n_comps) {
      tmp_comp_data <- comp_data %>% select(all_of(loadings %>% filter(component == j) %>% pull(feature)))
      tmp_inter_omic_cors <- get_all_inter_omic_corrs2(tmp_comp_data)
      mean_inter_omic_cors <- c(mean_inter_omic_cors, mean(tmp_inter_omic_cors$spearman_corr))
    }
    
    # Same, but using partial correlations to control for disease state
    mean_inter_omic_cors2 <- c()
    for (j in 1:n_comps) {
      tmp_comp_data <- comp_data %>% select(all_of(loadings %>% filter(component == j) %>% pull(feature)))
      tmp_inter_omic_cors <- get_all_inter_omic_corrs2(tmp_comp_data, control_for = data_Y_test)
      mean_inter_omic_cors2 <- c(mean_inter_omic_cors2, mean(tmp_inter_omic_cors$spearman_pcorr))
    }
    
    # Compute AUC for disease
    mean_disease_aucs <- c()
    for (j in 1:n_comps) {
      aucs <- c()
      for (l in 1:length(data_X)) {
        aucs <- c(aucs, as.numeric(auc(roc(quiet = T, response = data_Y_test, predictor = predictions[[l]][,j]))))
      }
      mean_disease_aucs <- c(mean_disease_aucs, mean(aucs))
    }
    
    # Record statistics
    component_associations <- bind_rows(
      component_associations,
      data.frame(
        method = 'DIABLO',
        design = round(diablo_design,2),
        perc_samples_excluded = perc_samples_to_exclude,
        n_samples_train = n_kept,
        keepX = diablo_keepX,
        i = i,
        comp_id = 1:n_comps,
        comp_size = c(table(loadings$component) %>% unname()),
        mean_inter_omic_cor = mean_inter_omic_cors,
        mean_inter_omic_pcor = mean_inter_omic_cors2,
        mean_disease_AUC = mean_disease_aucs
      )
    )
    
  } # Done iterating over random subsamples
  
} # Done iterating over keepX parameter

####################################################################
# MultiCCA
####################################################################

for (penalty_factor in c(1.8, 2)) { # penalty_factor <- 2 
  log_debug('penalty_factor = ', penalty_factor)
  
  # Run 100 times 
  for (i in 1:n_repeats) {
    
    # Randomly sample samples (rows)
    samples_subset <- sample(x = n_samples_tot, size = n_kept)
    samples_oof <- setdiff(1:n_samples_tot, samples_subset)
    data_X_train <- lapply(data_X, function(m, samps) return(m[samps,]), samps = samples_subset)
    data_Y_train <- data_Y[samples_subset]
    
    data_X_test <- lapply(data_X, function(m, samps) return(m[samps,]), samps = samples_oof)
    data_Y_test <- data_Y[samples_oof]
    
    # Run MultiCCA 
    tmp_cca_out1 <- MultiCCA(
      xlist = data_X_train,
      penalty = penalty_factor,
      niter = 50,
      ncomponents = n_comps,
      trace = FALSE
    )
    
    # Predict components for oof data
    predictions <- lapply(1:length(data_X_test), function(l) {
      tmp <- scale(data_X_test[[l]])
      tmp[is.na(tmp)] <- 0
      tmp %*% tmp_cca_out1$ws[[l]]
    })
    
    # Get loadings
    loadings <- get_loadings_smcca(data_X_train, tmp_cca_out1)
    
    # Compute correlations between variates for oof data
    comp_data <- bind_cols(lapply(data_X_test, as.data.frame))
    
    # Record mean correlations between features from different omics
    mean_inter_omic_cors <- c()
    for (j in 1:n_comps) {
      tmp_comp_data <- comp_data %>% select(all_of(loadings %>% filter(component == j) %>% pull(feature)))
      tmp_inter_omic_cors <- get_all_inter_omic_corrs2(tmp_comp_data)
      mean_inter_omic_cors <- c(mean_inter_omic_cors, mean(tmp_inter_omic_cors$spearman_corr))
    }
    
    # Same, but using partial correlations to control for disease state
    mean_inter_omic_cors2 <- c()
    for (j in 1:n_comps) {
      tmp_comp_data <- comp_data %>% select(all_of(loadings %>% filter(component == j) %>% pull(feature)))
      tmp_inter_omic_cors <- get_all_inter_omic_corrs2(tmp_comp_data, control_for = data_Y_test)
      mean_inter_omic_cors2 <- c(mean_inter_omic_cors2, mean(tmp_inter_omic_cors$spearman_pcorr))
    }
    
    # Compute AUC for disease
    mean_disease_aucs <- c()
    for (j in 1:n_comps) {
      aucs <- c()
      for (l in 1:length(data_X)) {
        aucs <- c(aucs, as.numeric(auc(roc(quiet = T, response = data_Y_test, predictor = predictions[[l]][,j]))))
      }
      mean_disease_aucs <- c(mean_disease_aucs, mean(aucs))
    }
    
    # Record statistics
    component_associations <- bind_rows(
      component_associations,
      data.frame(
        method = 'MultiCCA',
        perc_samples_excluded = perc_samples_to_exclude,
        n_samples_train = n_kept,
        keepX = penalty_factor,
        i = i,
        comp_id = 1:n_comps,
        comp_size = c(table(loadings$component) %>% unname()),
        mean_inter_omic_cor = mean_inter_omic_cors,
        mean_inter_omic_pcor = mean_inter_omic_cors2,
        mean_disease_AUC = mean_disease_aucs
      )
    )
    
  } # Done iterating over random subsamples
  
} # Done iterating over penalty parameter

####################################################################
# MintTea
####################################################################

for (diablo_keepX in c(10)) { # diablo_keepX <- 10 
  log_debug('diablo_keepX = ', diablo_keepX)
  
  # Run 100 times 
  for (i in 1:n_repeats) {
    set.seed(i)
    
    # Randomly sample samples (rows)
    samples_subset <- sample(x = n_samples_tot, size = n_kept)
    samples_oof <- setdiff(1:n_samples_tot, samples_subset)
    data_X_train <- lapply(data_X, function(m, samps) return(m[samps,]), samps = samples_subset)
    data_Y_train <- data_Y[samples_subset]
    
    data_X_test <- lapply(data_X, function(m, samps) return(m[samps,]), samps = samples_oof)
    data_Y_test <- data_Y[samples_oof]
    
    # Processed data
    proc_data_train <- data.frame(sample_id = samples_subset, disease_state = data_Y_train)
    for (ft in names(data_X_train)) proc_data_train <- bind_cols(proc_data_train, data_X_train[[ft]])
    
    # Run MintTea for data subset 1
    tmp_minttea_res <- MintTea(
      proc_data_train, 
      view_prefixes = names(data_X_train),
      param_diablo_keepX = c(diablo_keepX),
      param_diablo_design = c(diablo_design),
      param_edge_thresholds = c(minttea_edge_threshold)
      #return_main_results_only = FALSE
    )
    
    if (length(tmp_minttea_res) == 0) next;
    
    # Count number of final modules
    module_names <- names(tmp_minttea_res[[1]])
    n_modules <- length(module_names)

    # Compute correlations between variates for oof data
    proc_data_test <- bind_cols(lapply(data_X_test, as.data.frame))
    
    # Record mean correlations between features from different omics
    mean_inter_omic_cors <- c()
    for (j in module_names) {
      tmp_comp_data <- proc_data_test %>% select(all_of(tmp_minttea_res[[1]][[j]]$features))
      tmp_inter_omic_cors <- get_all_inter_omic_corrs2(tmp_comp_data)
      mean_inter_omic_cors <- c(mean_inter_omic_cors, mean(tmp_inter_omic_cors$spearman_corr))
    }
    
    # Same, but using partial correlations to control for disease state
    mean_inter_omic_cors2 <- c()
    for (j in module_names) {
      tmp_comp_data <- proc_data_test %>% select(all_of(tmp_minttea_res[[1]][[j]]$features))
      tmp_inter_omic_cors <- get_all_inter_omic_corrs2(tmp_comp_data, control_for = data_Y_test)
      mean_inter_omic_cors2 <- c(mean_inter_omic_cors2, mean(tmp_inter_omic_cors$spearman_pcorr))
    }
    
    # Predict variates for test data --> here we compute 1st PC's per view and per module
    latent_vars <- list()
    for (j in module_names) {
      pca_data <- proc_data_test %>% select(all_of(tmp_minttea_res[[1]][[j]]$features))
      pca_data <- pca_data[ , which(apply(pca_data, 2, var) != 0)] # Remove 0-variance columns
      if (ncol(pca_data) == 0) { latent_vars[[j]] <- NA; next }
      pca <- prcomp(pca_data, center = TRUE, scale. = TRUE) 
      latent_vars[[j]] <- pca$x[,'PC1']
    }
    
    # Compute AUC for disease
    aucs <- c()
    for (j in module_names) {
      if (all(is.na(latent_vars[[j]]))) next;
      aucs <- c(aucs, as.numeric(auc(roc(quiet = T, response = data_Y_test, predictor = latent_vars[[j]]))))
    }
    
    # Record statistics
    component_associations <- bind_rows(
      component_associations,
      data.frame(
        method = 'MintTea',
        design = round(diablo_design,2),
        perc_samples_excluded = perc_samples_to_exclude,
        n_samples_train = n_kept,
        keepX = diablo_keepX,
        i = i,
        comp_id = 1:n_modules,
        comp_size = sapply(tmp_minttea_res[[1]], function(x) { x$module_size }) %>% unname(),
        mean_inter_omic_cor = mean_inter_omic_cors,
        mean_inter_omic_pcor = mean_inter_omic_cors2,
        mean_disease_AUC = aucs
      )
    )
    
  } # Done iterating over random subsamples
  
} # Done iterating over keepX parameter

# View(component_associations %>% group_by(method, keepX) %>% summarize(mm = mean(comp_size)))

write_csv(component_associations, file = paste0("src/analyses/quality_comparisons_30perc_",d,".csv"))
message('DONE')
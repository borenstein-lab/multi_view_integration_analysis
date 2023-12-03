setwd('/specific/elhanan/PROJECTS/MULTI_VIEW_EM/repo/multi_view_integration_analysis')
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) stop("At least one argument must be supplied (dataset name).", call.=FALSE)
d <- args[1]
# d <- 'cd_franzosa_2019'
n_threads <- 32 # Threads

# Preperations
library(mixOmics)
library(dplyr)
library(logger)
library(PMA)
library(readr)
library(igraph)
set.seed(27)
options(scipen = 999)
log_threshold('DEBUG')

source("src/analyses/scca_and_diablo_utils.R")
source("src/intermediate_integration/MintTea.R")

# Load some dataset (after processing)
res_dir <- 'data/intermediate_integration_results'
f_path <- file.path(res_dir, paste0(d,'_results.RData')) 
load(f_path, envir = tmp_env <- new.env())
data_X <- tmp_env$minttea_results$diablo_input$X
data_Y <- tmp_env$minttea_results$diablo_input$Y
data_sample_ids <- tmp_env$minttea_results$diablo_input$sample_ids
n_samples_tot <- length(data_sample_ids)

# As components get less and less stable in CCA, we take only the first 3
n_comps <- 3

# Number of repeats
n_repeats <- 50

# Additional parameters
diablo_design <- .5 
minttea_edge_threshold <- .7
n_minttea_repeats <- 10

# Some utilities
get_latent_vars <- function(data_X, n_modules, modules_df) {
  latent_vars <- list()
  for (l in 1:length(data_X)) {
    latent_vars[[l]] <- list()
    for (j in 1:n_modules) {
      pca_data <- as.data.frame(data_X[[l]]) %>%
        select(any_of(modules_df %>% filter(module == paste0('module',j)) %>% pull(feature)))
      if(ncol(pca_data) == 0) { latent_vars[[l]][[j]] <- NA; next}
      pca <- prcomp(pca_data, center = TRUE, scale. = TRUE) 
      latent_vars[[l]][[j]] <- pca$x[,'PC1']
    }
  }
  return(latent_vars)
}

####################################################################
# DIABLO
####################################################################

# Record all identified modules for later debugging
debug <- data.frame()

# Results place holder
stability_results <- data.frame()

# Iterate over a few different diablo keepX params
for (diablo_keepX in c(7, 10)) { # diablo_keepX <- 7 
  log_debug('diablo_keepX = ', diablo_keepX)
  
  # Format to diablo's requirement
  list_keepX <- list()
  for(x in names(data_X)) list_keepX[[x]] <- rep(diablo_keepX, n_comps)
  
  # Iterate over percent samples to exclude
  for (perc_samples_to_exclude in c(0, 1, 5, 10, 20, 30, 50)) { # perc_samples_to_exclude <- 10 
    log_debug('perc_samples_to_exclude = ', perc_samples_to_exclude)
    
    n_excluded <- floor(n_samples_tot * perc_samples_to_exclude / 100)
    n_kept <- n_samples_tot - n_excluded
    
    # Run 100 times for each percent
    for (i in 1:n_repeats) {
      
      # Randomly sample once
      samples_subset_1 <- sample(x = n_samples_tot, size = n_kept)
      data_X_1 <- lapply(data_X, function(m, samps) return(m[samps,]), samps = samples_subset_1)
      data_Y_1 <- data_Y[samples_subset_1]
      # Shuffle order of columns too
      data_X_1 <- lapply(data_X_1, function(m) return(m[,sample(1:ncol(m))]))
      
      # Randomly sample twice
      samples_subset_2 <- sample(x = n_samples_tot, size = n_kept)
      data_X_2 <- lapply(data_X, function(m, samps) return(m[samps,]), samps = samples_subset_2)
      data_Y_2 <- data_Y[samples_subset_2]
      # Shuffle order of columns too
      data_X_2 <- lapply(data_X_2, function(m) return(m[,sample(1:ncol(m))]))
      
      # Run DIABLO for data subset 1
      tmp_diablo_out1 <- block.splsda(
        data_X_1, 
        data_Y_1, 
        ncomp = n_comps, 
        keepX = list_keepX, 
        max.iter = 300,
        design = diablo_design
      )
      
      # Extract components 1
      loadings_1 <- organize_diablo_loadings(tmp_diablo_out1) 
      
      # Run DIABLO for data subset 2
      tmp_diablo_out2 <- block.splsda(
        data_X_2, 
        data_Y_2, 
        ncomp = n_comps, 
        keepX = list_keepX, 
        max.iter = 300,
        design = diablo_design
      )
      
      # Extract components 2
      loadings_2 <- organize_diablo_loadings(tmp_diablo_out2) 
      
      # What % of features included in components 2 are also in component 1?
      n_intrsct <- n_distinct(intersect(loadings_2$feature,loadings_1$feature))
      perc_shared_feats <- 100 * n_intrsct / n_distinct(loadings_1$feature)
      jaccard <- n_intrsct / n_distinct(union(loadings_2$feature,loadings_1$feature))
      sorensen_dice <- 2 * n_intrsct / (n_distinct(loadings_1$feature) + n_distinct(loadings_2$feature))
        
      # Confusion matrix
      # all_feats <- sapply(data_X_1, colnames) %>% unlist(); matrix(c(n_distinct(intersect(loadings_2$feature,loadings_1$feature)), n_distinct(setdiff(loadings_1$feature, loadings_2$feature)), n_distinct(setdiff(loadings_2$feature, loadings_1$feature)), n_distinct(setdiff(x, union(loadings_1$feature, loadings_2$feature)))), nrow = 2)
      
      # What % of links included in components 2 are also in component 1? Calculate a z-score compared to shuffled components
      z_score_shared_links <- get_z_score_shared_links(loadings_1, loadings_2)
      
      # Are components orthogonal? if not - higher chances of more shared links
      mean_cors_within_comp1 <- c()
      for (l in 1:length(data_X)) {
        for (j in 2:n_comps) {
          for (k in 1:(j-1)) {
            mean_cors_within_comp1 <- c(mean_cors_within_comp1, abs(cor(tmp_diablo_out1$variates[[l]][,j], tmp_diablo_out1$variates[[l]][,k])))
          }
        }
      }
      mean_cors_within_comp1 <- mean(mean_cors_within_comp1)
      
      # Record statistics
      stability_results <- bind_rows(
        stability_results,
        data.frame(
          method = 'DIABLO',
          i = i,
          ncomp = n_comps,
          design = round(diablo_design,2),
          keepX = diablo_keepX,
          perc_samples_to_exclude = perc_samples_to_exclude,
          n_samples_kept = n_kept,
          n_feats_in_1 = nrow(loadings_1),
          n_feats_in_2 = nrow(loadings_2),
          mean_comp_size = table(loadings_1$component) %>% mean(),
          perc_shared_feats = perc_shared_feats,
          z_score_shared_links = z_score_shared_links,
          mean_cors_within_comp1 = mean_cors_within_comp1
        )
      )
      
      debug <- bind_rows(debug, loadings_1 %>% mutate(method='DIABLO',i=i,perc=perc_samples_to_exclude,penalty_factor=diablo_keepX,set='loadings1'))
      debug <- bind_rows(debug, loadings_2 %>% mutate(method='DIABLO',i=i,perc=perc_samples_to_exclude,penalty_factor=diablo_keepX,set='loadings2'))
      
    } # Done iterating over random subsamples
  
  } # Done iterating over percent to include

} # Done iterating over keepX parameter

####################################################################
# MultiCCA
####################################################################

# Iterate over a few different MultiCCA penalties
for (penalty_factor in c(1.8, 2)) { # penalty_factor <- 1 
  log_debug('penalty_factor = ', penalty_factor)
  
  # Iterate over percent samples to exclude
  for (perc_samples_to_exclude in c(0, 1, 5, 10, 20, 30, 50)) {  
    log_debug('perc_samples_to_exclude = ', perc_samples_to_exclude)
    
    n_excluded <- floor(n_samples_tot * perc_samples_to_exclude / 100)
    n_kept <- n_samples_tot - n_excluded
    
    # Run 100 times for each percent
    for (i in 1:n_repeats) {
      
      # Randomly sample once
      samples_subset_1 <- sample(x = n_samples_tot, size = n_kept)
      data_X_1 <- lapply(data_X, function(m, samps) return(m[samps,]), samps = samples_subset_1)
      # data_Y_1 <- data_Y[samples_subset_1]
      data_X_1 <- lapply(data_X_1, function(m) return(m[,sample(1:ncol(m))]))
      
      # Randomly sample twice
      samples_subset_2 <- sample(x = n_samples_tot, size = n_kept)
      data_X_2 <- lapply(data_X, function(m, samps) return(m[samps,]), samps = samples_subset_2)
      # data_Y_2 <- data_Y[samples_subset_2]
      data_X_2 <- lapply(data_X_2, function(m) return(m[,sample(1:ncol(m))]))
      
      # Run MultiCCA for data subset 1
      tmp_cca_out1 <- MultiCCA(
        xlist = data_X_1,
        penalty = penalty_factor,
        niter = 50,
        ncomponents = n_comps,
        trace = FALSE
      )
      
      # Extract components 1
      loadings_1 <- get_loadings_smcca(data_X_1, tmp_cca_out1)
      
      # Run MultiCCA for data subset 2
      tmp_cca_out2 <- MultiCCA(
        xlist = data_X_2,
        penalty = penalty_factor,
        niter = 50,
        ncomponents = n_comps,
        trace = FALSE
      )
      
      # Extract components 2
      loadings_2 <- get_loadings_smcca(data_X_2, tmp_cca_out2)
      
      # What % of features included in components 2 are also in component 1?
      n_intrsct <- n_distinct(intersect(loadings_2$feature,loadings_1$feature))
      perc_shared_feats <- 100 * n_intrsct / n_distinct(loadings_1$feature)
      jaccard <- n_intrsct / n_distinct(union(loadings_2$feature,loadings_1$feature))
      sorensen_dice <- 2 * n_intrsct / (n_distinct(loadings_1$feature) + n_distinct(loadings_2$feature))
      
      # What % of links included in components 2 are also in component 1? Calculate a z-score compared to shuffled components
      z_score_shared_links <- get_z_score_shared_links(loadings_1, loadings_2)
      
      # Are components orthogonal? if not - higher chances of more shared links
      mean_cors_within_comp1 <- c()
      for (l in 1:length(data_X)) {
        tmp1 <- scale(data_X_1[[l]]) %*% tmp_cca_out1$ws[[l]]
        for (j in 2:n_comps) {
          for (k in 1:(j-1)) {
            mean_cors_within_comp1 <- c(mean_cors_within_comp1, abs(cor(tmp1[,j], tmp1[,k])))
          }
        }
      }
      mean_cors_within_comp1 <- mean(mean_cors_within_comp1)
      
      # Record statistics
      stability_results <- bind_rows(
        stability_results,
        data.frame(
          method = 'MultiCCA',
          i = i,
          ncomp = n_comps,
          keepX = penalty_factor,
          perc_samples_to_exclude = perc_samples_to_exclude,
          n_samples_kept = n_kept,
          n_feats_in_1 = nrow(loadings_1),
          n_feats_in_2 = nrow(loadings_2),
          mean_comp_size = table(loadings_1$component) %>% mean(),
          perc_shared_feats = perc_shared_feats,
          z_score_shared_links = z_score_shared_links,
          mean_cors_within_comp1 = mean_cors_within_comp1
        )
      )
      
      debug <- bind_rows(debug, loadings_1 %>% mutate(method='MultiCCA',i=i,perc=perc_samples_to_exclude,penalty_factor=penalty_factor,set='loadings1'))
      debug <- bind_rows(debug, loadings_2 %>% mutate(method='MultiCCA',i=i,perc=perc_samples_to_exclude,penalty_factor=penalty_factor,set='loadings2'))
      
    } # Done iterating over random subsamples
    
  } # Done iterating over percent to include
  
} # Done iterating over penalty parameter

####################################################################
# MintTea
####################################################################

library(foreach)
library(doSNOW)
cluster <- makeCluster(n_threads) 
registerDoSNOW(cluster)

# For progress printing
pb <- txtProgressBar(max = n_repeats, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Iterate over a few different penalties
for (diablo_keepX in c(10)) { # diablo_keepX <- 10 
  log_debug('diablo_keepX = ', diablo_keepX)
  
  # Iterate over percent samples to exclude
  for (perc_samples_to_exclude in c(0, 1, 5, 10, 20, 30, 50)) {  
    log_debug('perc_samples_to_exclude = ', perc_samples_to_exclude)
    
    n_excluded <- floor(n_samples_tot * perc_samples_to_exclude / 100)
    n_kept <- n_samples_tot - n_excluded
    
    # Create a list of n_repeats*2 data subsamples
    samples_subsets_1 <- lapply(
      1:n_repeats, 
      function(i, n_samples_tot, n_kept) sample(x = n_samples_tot, size = n_kept), 
      n_samples_tot = n_samples_tot, 
      n_kept = n_kept)
    samples_subsets_2 <- lapply(
      1:n_repeats, 
      function(i, n_samples_tot, n_kept) sample(x = n_samples_tot, size = n_kept), 
      n_samples_tot = n_samples_tot, 
      n_kept = n_kept)
    
    # Run X times for each percent
    tmp_stability_results3 <- foreach(
      i = 1:n_repeats, 
      .combine = rbind, 
      .packages = c('mixOmics', 'igraph', 'dplyr'), 
      .options.snow = opts) %dopar% {
      
      # Randomly sample once
      samples_subset_1 <- samples_subsets_1[[i]]
      data_X_1 <- lapply(data_X, function(m, samps) return(m[samps,]), samps = samples_subset_1)
      data_Y_1 <- data_Y[samples_subset_1]
      data_X_1 <- lapply(data_X_1, function(m) return(m[,sample(1:ncol(m))])) # Shuffle columns
      
      # Randomly sample twice
      samples_subset_2 <- samples_subsets_2[[i]]
      data_X_2 <- lapply(data_X, function(m, samps) return(m[samps,]), samps = samples_subset_2)
      data_Y_2 <- data_Y[samples_subset_2]
      data_X_2 <- lapply(data_X_2, function(m) return(m[,sample(1:ncol(m))])) # Shuffle columns
      
      # Run MintTea for data subset 1
      proc_data_1 <- data.frame(sample_id = samples_subset_1, disease_state = data_Y_1)
      for (ft in names(data_X_1)) proc_data_1 <- bind_cols(proc_data_1, data_X_1[[ft]])
      tmp_minttea_res_1 <- MintTea(
        proc_data_1, 
        view_prefixes = names(data_X_1),
        param_diablo_keepX = c(diablo_keepX),
        param_diablo_design = c(diablo_design),
        param_edge_thresholds = c(minttea_edge_threshold)
      )
      
      # Run MintTea for data subset 2
      proc_data_2 <- data.frame(sample_id = samples_subset_2, disease_state = data_Y_2)
      for (ft in names(data_X_2)) proc_data_2 <- bind_cols(proc_data_2, data_X_2[[ft]])
      tmp_minttea_res_2 <- MintTea(
        proc_data_2, 
        view_prefixes = names(data_X_2),
        param_diablo_keepX = c(diablo_keepX),
        param_diablo_design = c(diablo_design),
        param_edge_thresholds = c(minttea_edge_threshold)
      )
      
      # Comparison irrelavant if no modules detected at all
      if (length(tmp_minttea_res_1) == 0 | length(tmp_minttea_res_2) == 0) {
        n_modules <- ifelse(length(tmp_minttea_res_1) == 0, 0, length(tmp_minttea_res_1[[1]]))
        data.frame(                                       
          method = 'MintTea',
          i = i,
          ncomp = n_modules,
          keepX = diablo_keepX,
          perc_samples_to_exclude = perc_samples_to_exclude,
          n_samples_kept = n_kept,
          n_feats_in_1 = NA,
          n_feats_in_2 = NA,
          mean_comp_size = NA,
          perc_shared_feats = NA,
          z_score_shared_links = NA,
          mean_cors_within_comp1 = NA # mean_cors_within_comp1
        )
      } else {
        # What % of features included in components 2 are also in component 1?
        all_feats_1 <- sapply(tmp_minttea_res_1[[1]], function(x) x$features) %>% unlist() %>% unname() 
        all_feats_2 <- sapply(tmp_minttea_res_2[[1]], function(x) x$features) %>% unlist() %>% unname() 
        n_intrsct <- length(intersect(all_feats_1,all_feats_2))
        perc_shared_feats <- 100 * n_intrsct / length(all_feats_1)
        jaccard <- n_intrsct / n_distinct(union(all_feats_1, all_feats_2))
        sorensen_dice <- 2 * n_intrsct / (length(all_feats_1) + length(all_feats_2))
        
        # What % of links included in components 2 are also in component 1? Calculate a z-score compared to shuffled components
        module_names_1 <- names(tmp_minttea_res_1[[1]])
        modules_1 <- bind_rows(lapply(
          module_names_1, 
          function(x, tmp_minttea_res) { return(data.frame(feature = tmp_minttea_res[[1]][[x]]$features, module = x)) },
          tmp_minttea_res = tmp_minttea_res_1
        )) %>% mutate(feature_set = substr(feature,0,1))
        module_names_2 <- names(tmp_minttea_res_2[[1]])
        modules_2 <- bind_rows(lapply(
          module_names_2, 
          function(x, tmp_minttea_res) { return(data.frame(feature = tmp_minttea_res[[1]][[x]]$features, module = x)) },
          tmp_minttea_res = tmp_minttea_res_2
        )) %>% mutate(feature_set = substr(feature,0,1))
        z_score_shared_links <- get_z_score_shared_links(modules_1, modules_2, join_by = 'module')
        
        
        # Count number of final modules
        n_modules <- n_distinct(modules_1$module)
        
        # # Predict variates for train data --> 1st PC's per view and per module
        # latent_vars1 <- get_latent_vars(data_X_1, n_modules, modules_1)
        # latent_vars2 <- get_latent_vars(data_X_2, n_modules, modules_2)
        # 
        # # Are components orthogonal? if not - higher chances of more shared links
        # if (n_modules == 1) {
        #   mean_cors_within_comp1 <- NA
        #   } else {
        #   mean_cors_within_comp1 <- c()
        #   for (l in 1:length(data_X)) {
        #     for (j in 2:n_modules) {
        #       for (k in 1:(j-1)) {
        #         if ((! all(is.na(latent_vars1[[l]][[j]]))) & (! all(is.na(latent_vars1[[l]][[k]]))))
        #           mean_cors_within_comp1 <- c(mean_cors_within_comp1, abs(cor(latent_vars1[[l]][[j]], latent_vars1[[l]][[k]])))
        #       }
        #     }
        #   }
        #   mean_cors_within_comp1 <- mean(mean_cors_within_comp1)
        # }
        
        # Record statistics
        data.frame(                                       
          method = 'MintTea',
          i = i,
          ncomp = n_modules,
          keepX = diablo_keepX,
          perc_samples_to_exclude = perc_samples_to_exclude,
          n_samples_kept = n_kept,
          n_feats_in_1 = nrow(modules_1),
          n_feats_in_2 = nrow(modules_2),
          mean_comp_size = table(modules_1$module) %>% mean(),
          perc_shared_feats = perc_shared_feats,
          z_score_shared_links = z_score_shared_links,
          mean_cors_within_comp1 = NA # mean_cors_within_comp1
        )
      } # End of "else"
      
    } # Done iterating over random subsamples
    
    stability_results <- bind_rows(stability_results, tmp_stability_results3)
    
  } # Done iterating over percent to include
  
} # Done iterating over penalty parameter

stopCluster(cluster)
write_csv(stability_results, file = paste0("src/analyses/stability_analysis_",d,".csv"))
message('DONE')

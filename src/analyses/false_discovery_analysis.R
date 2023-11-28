setwd('/specific/elhanan/PROJECTS/MULTI_VIEW_EM/repo/multi_view_integration_analysis')
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) stop("At least one argument must be supplied (dataset name).", call.=FALSE)
d <- args[1]
# d <- 'cd_franzosa_2019'

# Preperations
library(mixOmics)
library(dplyr)
library(logger)
library(PMA)
library(readr)
library(igraph)
set.seed(2727)
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

# As components get less and less stable, we take only the first 3
n_comps <- 3

# Number of repeats
n_repeats <- 100

# Two shuffling approaches
shuffle_within_omic <- F

# We use a diablo_design as used in the stability analysis
diablo_design <- .5
minttea_edge_threshold <- .7

# Results place holder
false_discovery_results <- data.frame()

####################################################################
# DIABLO
####################################################################

# Iterate over a few different diablo keepX params
for (diablo_keepX in c(7,10)) { # diablo_keepX <- 10 
  log_debug('DIABLO: diablo_keepX = ', diablo_keepX)
  
  # Format to diablo's requirement
  list_keepX <- list()
  for(x in names(data_X)) list_keepX[[x]] <- rep(diablo_keepX, n_comps)
  
  # Run n_repeats times 
  for (i in 1:n_repeats) {
    
    # Shuffle data to "break" inter-omic associations
    data_X_shuf <- lapply(data_X, function(m) return(m[sample(1:nrow(m)),]))
    
    if (shuffle_within_omic) { data_X_shuf = lapply(data_X_shuf, function(m) return(apply(m, 2, sample))) }
    
    # Run DIABLO for shuffled data
    tmp_diablo_out <- block.splsda(
      data_X_shuf, 
      data_Y, 
      ncomp = n_comps, 
      keepX = list_keepX, 
      max.iter = 300,
      design = diablo_design
    )
    
    # Evaluate how well latent variables correlate and discard non-significantly correlated components
    signif_comps <- c()
    for (j in 1:n_comps) {
      ps <- c()
      for (l in 1:(length(data_X)-1)) {
        for (k in (l+1):length(data_X)) {
          tmp_cor <- cor.test(tmp_diablo_out$variates[[l]][,j], tmp_diablo_out$variates[[k]][,j])
          ps <- c(ps, tmp_cor$p.value)
        }
      }
      ps <- p.adjust(ps, method = 'bonf')
      if(all(ps < 0.01)) signif_comps <- c(signif_comps, j)
    }
    
    # Get loadings
    loadings <- organize_diablo_loadings(tmp_diablo_out) %>%
      filter(component %in% signif_comps) %>%
      filter(abs(loading) > 0.001)
    
    # How many modules (multi-omic?) did we find? Record statistics
    false_discovery_results <- bind_rows(
      false_discovery_results,
      data.frame(
        method = 'DIABLO',
        i=i,
        original_ncomp = n_comps,
        design = round(diablo_design,2),
        keepX = diablo_keepX,
        final_ncomp = length(signif_comps),
        mean_comp_size = ifelse(length(signif_comps) == 0, 0, table(loadings$component) %>% mean() %>% round(2))
      )
    )
    
  } # Done iterating over random subsamples
  
} # Done iterating over keepX parameter

####################################################################
# MultiCCA
####################################################################

# Iterate over a few different params
for (penalty_factor in c(1.8, 2)) { 
  log_debug('MultiCCA: penalty_factor = ', penalty_factor)
  
  # Format to diablo's requirement
  list_keepX <- list()
  for(x in names(data_X)) list_keepX[[x]] <- rep(diablo_keepX, n_comps)
  
  # Run n_repeats times 
  for (i in 1:n_repeats) {
    
    # Shuffle data to "break" inter-omic associations
    data_X_shuf <- lapply(data_X, function(m) return(m[sample(1:nrow(m)),]))
    
    if (shuffle_within_omic) { data_X_shuf = lapply(data_X_shuf, function(m) return(apply(m, 2, sample))) }
    
    tmp_cca_out1 <- MultiCCA(
      xlist = data_X_shuf,
      penalty = penalty_factor,
      niter = 50,
      ncomponents = n_comps,
      trace = FALSE
    )
    
    # Evaluate how well latent variables correlate and discard non-significantly correlated components
    signif_comps <- c()
    for (j in 1:n_comps) {
      ps <- c()
      for (l in 1:(length(data_X)-1)) {
        for (k in (l+1):length(data_X)) {
          tmp_cor <- cor.test(scale(data_X_shuf[[l]]) %*% tmp_cca_out1$ws[[l]][,j], scale(data_X_shuf[[k]]) %*% tmp_cca_out1$ws[[k]][,j])
          ps <- c(ps, tmp_cor$p.value)
        }
      }
      ps <- p.adjust(ps, method = 'bonf')
      if(all(ps < 0.01)) signif_comps <- c(signif_comps, j)
    }
    
    # Get loadings
    loadings <- get_loadings_smcca(data_X_shuf, tmp_cca_out1) %>%
      filter(component %in% signif_comps) %>%
      filter(abs(ws) > 0.001)
    
    # How many modules (multi-omic?) did we find? Record statistics
    false_discovery_results <- bind_rows(
      false_discovery_results,
      data.frame(
        method = 'MultiCCA',
        i=i,
        original_ncomp = n_comps,
        design = round(diablo_design,2),
        keepX = penalty_factor,
        final_ncomp = length(signif_comps),
        mean_comp_size = ifelse(length(signif_comps) == 0, 0, table(loadings$component) %>% mean() %>% round(2))
      )
    )
    
  } # Done iterating over random subsamples
} # Done iterating over keepX parameter

####################################################################
# MintTea
####################################################################

# Iterate over a few different diablo keepX params
for (diablo_keepX in c(10)) { # diablo_keepX <- 10 
  log_debug('MintTea: diablo_keepX = ', diablo_keepX)
  
  # Format to diablo's requirement
  list_keepX <- list()
  for(x in names(data_X)) list_keepX[[x]] <- rep(diablo_keepX, n_comps)
  
  # Run n_repeats times 
  for (i in 1:n_repeats) {
    set.seed(i)
    
    # Shuffle data to "break" inter-omic associations
    data_X_shuf <- lapply(data_X, function(m) return(m[sample(1:nrow(m)),]))
    
    if (shuffle_within_omic) { data_X_shuf = lapply(data_X_shuf, function(m) return(apply(m, 2, sample))) }
    
    # Processed data
    proc_data <- data.frame(sample_id = data_sample_ids, disease_state = data_Y)
    for (ft in names(data_X_shuf)) proc_data <- bind_cols(proc_data, data_X_shuf[[ft]])
    
    # Run MintTea for data subset 1
    tmp_minttea_res <- MintTea(
      proc_data, 
      view_prefixes = names(data_X_shuf),
      param_diablo_keepX = c(diablo_keepX),
      param_diablo_design = c(diablo_design),
      param_edge_thresholds = c(minttea_edge_threshold)
      #return_main_results_only = FALSE
    )
    
    # How many modules (multi-omic?) did we find? Record statistics
    if (length(tmp_minttea_res) == 0) {
      false_discovery_results <- bind_rows(
        false_discovery_results,
        data.frame(
          method = 'MintTea',
          i=i,
          original_ncomp = n_comps,
          design = round(diablo_design,2),
          keepX = diablo_keepX,
          final_ncomp = 0,
          mean_comp_size = 0
        )
      )
    } else {
      false_discovery_results <- bind_rows(
        false_discovery_results,
        data.frame(
          method = 'MintTea',
          i=i,
          original_ncomp = n_comps,
          design = round(diablo_design,2),
          keepX = diablo_keepX,
          final_ncomp = length(tmp_minttea_res[[1]]),
          mean_comp_size = sapply(tmp_minttea_res[[1]], function(x) { x$module_size }) %>% mean() %>% round(2)
        )
      )
    }

  } # Done iterating over random subsamples
} # Done iterating over keepX parameter

write_csv(false_discovery_results, file = paste0("src/analyses/false_discovery_results_",d,".csv"))
message('DONE')

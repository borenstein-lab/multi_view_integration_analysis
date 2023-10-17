########################################################################
# This script contains the main intermediate integration pipeline used 
#  for the manuscript, using MintTea (Multi-view INTegration Tool for 
#  MicrobiomE Analysis).
# 
# To run MintTea on your own data, please refer to MintTea.R script.
########################################################################

########################################################################
# Read command line arguments
########################################################################

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

########################################################################
# Preparations
########################################################################

library(config)
library(dplyr)
library(readr)

source('src/ml_pipeline/clustering.R')
source("src/ml_pipeline/preprocessing.R")
source("src/intermediate_integration/MintTea.R")

# Read configuration file
config <- config::get(file = "src/ml_pipeline/config.yml")

log_threshold('DEBUG')
options(scipen = 999)

########################################################################
# Define output files 
########################################################################

res_dir <- 'data/intermediate_integration_results'
dir.create(res_dir, showWarnings = FALSE)
out_clusters <- paste0(res_dir, '/', d, "_clusters.csv") # Ignore
out_rdata <- paste0(res_dir, '/', d, '_results.RData')

########################################################################
# Load and organize data 
########################################################################

# Load processed data
log_debug("Loading dataset: {d}")
include_metabs <- file.exists(sprintf(config$paths_templates$metabolites, d))

# File names
metadata_path <- sprintf(config$paths_templates$metadata, d)
taxonomy_path <- sprintf(config$paths_templates$taxonomy, d)
pathways_path <- sprintf(config$paths_templates$pathways, d)
if (include_metabs) metabolites_path <- sprintf(config$paths_templates$metabolites, d)
metagenome_path <- sprintf(config$paths_templates$metagenome, d)

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

########################################################################
# Intermediate integration using MintTea
########################################################################

minttea_results <- MintTea(
  proc_data, 
  view_prefixes = c('T', 'G', 'P', 'M'), 
  return_main_results_only = FALSE
)

# Organize and save results in an RData file
save(minttea_results,file = out_rdata)
log_debug('Saved MintTea objects to: {out_rdata}')


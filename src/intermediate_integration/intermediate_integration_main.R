########################################################################
# This script contains the main intermediate integration pipeline used 
#  for the manuscript.
# 
# To run MintTea on your own data, please refer to MintTea.R script.
########################################################################

# udocker run --volume=/specific/elhanan/PROJECTS/MULTI_VIEW_EM/repo/multi_view_integration_analysis efrat_ubun_r Rscript /specific/elhanan/PROJECTS/MULTI_VIEW_EM/repo/multi_view_integration_analysis/src/intermediate_integration/intermediate_integration_main.R

# setwd('/specific/elhanan/PROJECTS/MULTI_VIEW_EM/repo/multi_view_integration_analysis')


# Preparations
library(dplyr)
library(readr)

source("src/intermediate_integration/MintTea.R")

log_threshold('DEBUG')
options(scipen = 999)

# Output folders
res_dir <- 'data/intermediate_integration_results'
dir.create(res_dir, showWarnings = FALSE)

# Datasets and views/omics
DATASETS <- list(
  'crc_zeller_2014' = c('T','P'),
  'crc_vogtmann_2016' = c('T','P'),
  'crc_wirbel_2018' = c('T','P'),
  'igt_karlsson_2013' = c('T','P'),
  't2d_karlsson_2013' = c('T','P'),
  'schizofrenia_zhu_2020' = c('T','P'),
  'me_cfs_nagySzakal_2017' = c('T','P'),
  'ht_li_2017' = c('T','P'),
  'pre_ht_li_2017' = c('T','P'),
  # ---
  'cirrhosis_qin_2014' = c('T','P'),
  'crc_feng_2015' = c('T','P'),
  'crc_yu_2015' = c('T','P'),
  'sth_rubel_2020' = c('T','P'),
  'uc_spain_nielsen_2014' = c('T','P'),
  'cd_franzosa_2019' = c('T','P','M'),
  'crc_s3_s4_yachida_2019' = c('T','P','M'),
  'esrd_wang_2020' = c('T','P','M'),
  'uc_franzosa_2019' = c('T','P','M'),
  'metacardis_1_8' = c('T','P','S'),
  'metacardis_3_8' = c('T','P','S')
)

for (d in names(DATASETS)) {
  
  # Set file names
  input_data <- file.path('data/ml_input',d, 'all_data.tsv')
  out_rdata <- paste0(res_dir, '/', d, '_results.RData')
  
  # Load processed data
  log_debug("Loading dataset: {d}")
  proc_data <- read_delim(input_data, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
  
  # Rename for simplicity
  if ('DiseaseState' %in% names(proc_data)) proc_data <- proc_data %>% rename(disease_state = DiseaseState)
  if ('sample_id__' %in% names(proc_data)) proc_data <- proc_data %>% rename(sample_id = sample_id__)

  # Intermediate integration using MintTea
  minttea_results <- MintTea(
    proc_data, 
    view_prefixes = DATASETS[[d]], 
    return_main_results_only = FALSE,
    param_diablo_keepX = c(10,15),
    param_diablo_design = c(0.3,0.5),
    param_edge_thresholds = c(0.7,0.8)
  )

  # View(minttea_results$modules_list)

  # Organize and save results in an RData file
  save(minttea_results,file = out_rdata)
  log_debug('Saved MintTea objects to: {out_rdata}')
}

log_debug('Done')

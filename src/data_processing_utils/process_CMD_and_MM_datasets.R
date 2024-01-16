library(readr)
library(tidyr)
library(dplyr)
library(logger)

source('src/ml_pipeline/clustering.R')
source('src/ml_pipeline/preprocessing.R')
log_threshold(DEBUG)

CMD_DATASETS <- c(
  'cirrhosis_qin_2014',
  'crc_feng_2015',
  'crc_yu_2015',
  'sth_rubel_2020',
  'uc_spain_nielsen_2014',
  # ---
  'crc_zeller_2014',
  'crc_vogtmann_2016',
  'crc_wirbel_2018',
  'igt_karlsson_2013',
  't2d_karlsson_2013',
  'schizofrenia_zhu_2020',
  'me_cfs_nagySzakal_2017',
  'ht_li_2017',
  'pre_ht_li_2017'
)

MM_DATASETS <- c(
  'cd_franzosa_2019',
  'crc_s3_s4_yachida_2019',
  'esrd_wang_2020',
  'uc_franzosa_2019'
)

for (d in c(CMD_DATASETS, MM_DATASETS)) { # d = 'cd_franzosa_2019'
  log_debug('Preparing dataset: ', d)
  
  # Paths for input files
  pwy_file <- file.path('data/ml_input',d,'pathways.tsv')
  tax_file <- file.path('data/ml_input',d,'taxonomy.tsv')
  metadata_file <- file.path('data/ml_input',d,'metadata.tsv')
  if (d %in% MM_DATASETS) mtb_file <- file.path('data/ml_input',d,'metabolites.tsv')
    
  # Paths for output files
  clusters_file <- file.path('data/ml_input',d,'clusters.tsv')
  proc_data_file <- file.path('data/ml_input',d,'all_data.tsv')
  
  # -----------------------------------------------
  # Organize metadata 
  # -----------------------------------------------
  
  metadata <- read_delim(metadata_file, delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE, show_col_types = F) %>%
    rename(sample_id__ = sample_id) %>%
    mutate(sample_id__ = as.character(sample_id__)) %>%
    mutate(DiseaseState = ifelse(DiseaseState == 'H', 'healthy', 'disease'))

  # -----------------------------------------------
  # Organize species table (add T__ prefix)
  # -----------------------------------------------
  
  species <- prep_load_taxonomy(tax_file) 

  # Remove rare features + constant features
  species <- prep_sanitize_dataset(species, 'T', rare_feature_cutoff = 0.15, mean_abundance_cutoff = 0.00005)

  # Some correlation statistics
  # ggplot(species, aes(x = T__Oscillibacter.ruminantium, y = T__Oscillibacter.valericigenes)) + geom_point(alpha = 0.3) + theme_bw() + ggtitle('Correlation between relative abundances') + theme(plot.title = element_text(hjust = 0.5))
  
  # -----------------------------------------------
  # Organize pathways table (add P__ prefix)
  # -----------------------------------------------

  pwy <- prep_load_pathways(pwy_file)

  # Remove rare features + constant features
  pwy <- prep_sanitize_dataset(pwy, 'P', rare_feature_cutoff = 0.15, mean_abundance_cutoff = 0.00005)

  # -----------------------------------------------
  # Organize metabolomics (add M__ prefix)
  # -----------------------------------------------

  if (d %in% MM_DATASETS) {
    mtb <- prep_load_metabolites(mtb_file)
    
    # Remove rare features + constant features
    mtb <- prep_sanitize_dataset(mtb %>% mutate_if(is.numeric, function(x) ifelse(is.na(x), 0, x)), 'M')
    
    # Log transform metabolite values
    # (First - replace zeros/NAs with half the minimal non-zero value)
    mtb <- mtb %>%
      mutate_if(is.numeric, function(x) ifelse(x==0, NA, x)) %>%
      mutate_if(is.numeric, function(x) ifelse(is.na(x), min(x, na.rm = T)/2, x)) %>%
      mutate_if(is.numeric, log)
  }
  
  # -----------------------------------------------
  # Cluster highly correlated features
  # -----------------------------------------------

  all_data <- metadata %>%
    select(sample_id__, DiseaseState) %>%
    inner_join(species, by = 'sample_id__') %>%
    inner_join(pwy, by = 'sample_id__') 
    
  if (d %in% MM_DATASETS) 
    all_data <- all_data %>% inner_join(mtb, by = 'sample_id__')

  cluster_init_clusters_table(clusters_file)

  all_data_clustered <- cluster_cluster_features(
    dataset = all_data %>% select(-sample_id__),
    feature_set_type = 'All',
    cluster_type = 'clustering99',
    clusters_output = clusters_file
  )
  all_data <- bind_cols(all_data %>% select(sample_id__), all_data_clustered)

  # ---------------------------------------------------
  # Save to files
  # ---------------------------------------------------

  # Save for intermediate integration (one large table)
  write_tsv(all_data, proc_data_file)

  # Print summaries (numbers of cases/controls)
  log_debug('Labels summary: ', 
            paste(names(table(all_data$DiseaseState)), collapse = ','), ' -- ', 
            paste(table(all_data$DiseaseState), collapse = ','))
  
}

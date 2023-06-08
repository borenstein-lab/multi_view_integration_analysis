library(config)
library(logger)
library(dplyr)
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

source('pathways_db.R')

prep_add_col_prefix <- function(df, prefix) {
  sample_id_col <- df %>% select(sample_id__)
  df <- df %>% select(-sample_id__)
  colnames(df) <- paste0(prefix, "__", colnames(df))
  df <- bind_cols(sample_id_col, df)
  return(df)
}

# Performs TSS normalization
prep_convert_to_relative_abundances <- function(df) {
  sample_id_col <- df %>%
    select(sample_id__) %>%
    mutate(sample_id__ = as.character(sample_id__))
  
  df <- df %>%
    select(-sample_id__) %>%
    mutate_all(as.numeric) %>%
    mutate(totalSum = rowSums(across(everything()))) %>%
    mutate(across(.cols = -c(totalSum), ~ .x / totalSum)) %>%
    select(-totalSum) %>%
    bind_cols(sample_id_col)
  return(df)
}

# Performs a CLR transformation on a given data table
prep_clr_transpose <- function(df, pseudo_count) {
  require(mixOmics)
  includes_sample_id <- ifelse('sample_id__' %in% colnames(df), TRUE, FALSE)
  tmp_df <- df
  if(includes_sample_id) tmp_df <- df %>% select(-sample_id__)
  
  # CLR transform values
  clr_df <- logratio.transfo(tmp_df, logratio = 'CLR', offset = pseudo_count)
  
  # Reformat as table and add sample ids if needed
  class(clr_df) <- "matrix"
  if (!includes_sample_id) return(data.frame(clr_df))
  return(bind_cols(data.frame(clr_df), df %>% select(sample_id__)))
}

prep_load_metadata <- function(path) {
  valid_labels <- c('H', config::get('disease_labels'))

  metadata_df <- read.csv(file = path, sep = '\t', header = TRUE) %>%
      as_tibble() %>%
      rename(sample_id__ = 1) %>%
      mutate(sample_id__ = as.character(sample_id__)) %>%
      select(sample_id__, DiseaseState) %>%
      filter(DiseaseState %in% valid_labels) %>%
      mutate(DiseaseState = replace(DiseaseState, DiseaseState == 'H', 'healthy')) %>%
      mutate(DiseaseState = replace(DiseaseState, DiseaseState != 'healthy', 'disease')) %>%
      mutate_at(vars(DiseaseState), factor)
  return(metadata_df)
}


prep_load_pathways <- function(path, 
                               remove_non_bacterial = TRUE, 
                               remove_super_pwys = FALSE, 
                               convert_to_relab = TRUE) {
  # Load pathways and filter non-bacterial pathways and redundant (super-pathways) pathways
  pathways_df <- read.table(file = path, sep = '\t',
                            header = TRUE,
                            quote = "",
                            check.names = FALSE) 
    
  if (remove_non_bacterial) pathways_df <- filter_bacterial_pathways(pathways_df)
  if (remove_super_pwys) pathways_df <- filter_super_pathways(pathways_df)
  
  # Transpose pathways DF and correct column names to pathways names
  pathways_df <- pathways_df %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column() %>%
    as_tibble()

  names(pathways_df) <- unlist(pathways_df[1, ])
  
  pathways_df <- pathways_df[-(1:2), ]
  pathways_df <- rename(pathways_df, sample_id__ = 1)
  pathways_df <- prep_add_col_prefix(pathways_df, 'P')
  
  if (convert_to_relab) pathways_df <- prep_convert_to_relative_abundances(pathways_df)
  return(pathways_df)
}

prep_load_metagenome <- function(path, convert_to_relab = TRUE) {
  # Read table, remove "description" column if exists (outputted by picrust)
  df <- read.csv(file = path, header = TRUE, check.names = FALSE, sep = '\t') %>%
    select(-any_of(c('description'))) 
  
  if ("function" %in% names(df)) 
    df <- df %>% rename(`KO` = `function`)
  
  df <- df %>%
    tibble::column_to_rownames(var = "KO") %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "sample_id__") %>%
    prep_add_col_prefix('G')
  
  if (convert_to_relab) df <- prep_convert_to_relative_abundances(df)
  return(df)
}

prep_load_taxonomy <- function(path, convert_to_relab = TRUE) {
  df <- read.csv(file = path, sep = '\t', header = TRUE, check.names = FALSE) %>%
    select(-any_of(config::get('unknown_taxa'))) %>%
    rename(sample_id__ = 1) %>%
    prep_add_col_prefix('T')
  
  if (convert_to_relab) df <- prep_convert_to_relative_abundances(df)
  return(df)
}


prep_load_metabolites <- function(path, remove_unknown_metabolites = TRUE) {
  metabolites_df <- read.csv(file = path, sep = '\t', header = TRUE, check.names = FALSE) %>%
    rename(sample_id__ = 1) %>%
    mutate(sample_id__ = as.character(sample_id__)) %>%  
    prep_add_col_prefix('M')
  
  # For specific datasets with true metabolomics - we only keep identified metabolites to facilitate interpretation
  d <- basename(dirname(path))
  if (remove_unknown_metabolites & d == "esrd_wang_2020") {
    keep <- grep("_unknown$", colnames(metabolites_df), invert = TRUE, value = TRUE)
    metabolites_df <- metabolites_df[,keep]
  } else if (remove_unknown_metabolites & d %in% c("uc_franzosa_2019","cd_franzosa_2019")) {
    keep <- grep("\\: NA$", colnames(metabolites_df), invert = TRUE, value = TRUE) 
    metabolites_df <- metabolites_df[,keep]
  }
  
  return(metabolites_df)
}


prep_join_metadata <- function(dataset_df, metadata_df) {
  all_data <- dataset_df %>%
    inner_join(metadata_df, by = 'sample_id__')
  return(all_data)
}


prep_join_features <- function(df1, df2) {
  joined_df <- df1 %>% 
    inner_join(df2, by = 'sample_id__') %>%
    relocate(sample_id__)
  return(joined_df)
}


prep_sanitize_dataset <- function(df, feature_set, rare_feature_cutoff = 0.15, low_abundance_cutoff = 0.00005) {
  
  # Remove constant features
  nfeatures <- ncol(df)
  df <- df %>% select(where(~ n_distinct(.) > 1))
  log_debug(sprintf('Sanitizer: removed %d/%d (%d%%) constant features for dataset %s',
                    (nfeatures - ncol(df)),
                    nfeatures,
                    as.integer(100 * (nfeatures - ncol(df)) / nfeatures),
                    feature_set))
  
  # Remove rare features (have <15% non-zero values)
  nfeatures <- ncol(df)
  non_zero_percentage <- colSums(df != 0) / nrow(df)
  rare_features <- names(non_zero_percentage[non_zero_percentage <= rare_feature_cutoff])
  df <- df %>% select(-all_of(rare_features))
  log_debug(sprintf('Sanitizer: removed %d/%d (%d%%) low-prevalance features for dataset %s',
                    (nfeatures - ncol(df)),
                    nfeatures,
                    as.integer(100 * (nfeatures - ncol(df)) / nfeatures),
                    feature_set))
  
  # Remove rare features by mean abundance (only for TSS-normalized feature types)
  if (feature_set %in% c('T','G','P')) {
    # Verify the data looks like TSS-normalized data, i.e. each row sums to ~1.
    if (sum(round(rowSums(df %>% select(-sample_id__))) == 1) / nrow(df) < 0.95)
      log_error('Seems like the feature table is not TSS-normalized. Abundance filtering should not be performed')
    nfeatures <- ncol(df)
    mean_abundances <- apply(df %>% select(-sample_id__), 2, mean)
    rare_features <- names(mean_abundances[mean_abundances <= low_abundance_cutoff])
    df <- df %>% select(-all_of(rare_features))
    log_debug(sprintf('Sanitizer: removed %d/%d (%d%%) low-abundance features for dataset %s',
                      (nfeatures - ncol(df)),
                      nfeatures,
                      as.integer(100 * (nfeatures - ncol(df)) / nfeatures),
                      feature_set))
  }
  
  return(df)
}


# There is an option for override configuration clustering, used by cross.
# If NULL is supplied, the clustering config is used. Otherwise, the clustering
# specifid by "override_clustering" arg is used.
prep_preprocess_dataset <- function(ds_name, 
                                    cluster_type = config::get('cluster_type'),
                                    feature_set_types = c('T','G','P','M','T+G+P','T+G+P+M'), 
                                    clr_trans = FALSE) {
  # Get file paths
  metadata_path <- sprintf(config::get('paths_templates')$metadata, ds_name)
  taxonomy_path <- sprintf(config::get('paths_templates')$taxonomy, ds_name)
  pathways_path <- sprintf(config::get('paths_templates')$pathways, ds_name)
  metabolites_path <- sprintf(config::get('paths_templates')$metabolites, ds_name)
  genes_path <- sprintf(config::get('paths_templates')$metagenome, ds_name)
  clusters_output <- sprintf(config::get('paths_templates')$clusters, ds_name)
  
  # If there are no metabolomics available, update list of feature sets
  include_metabs <- ifelse(file.exists(metabolites_path), TRUE, FALSE)
  if(!include_metabs) {
    feature_set_types <- feature_set_types[!grepl('M', feature_set_types)]
    log_debug("No metabolomics in this dataset - removing M feature set type")
  } else {
    feature_set_types <- feature_set_types[feature_set_types != 'T+G+P']
  }
  
  # Initialize clusters file
  cluster_init_clusters_table(clusters_output)

  # Load datasets
  metadata_df <- prep_load_metadata(metadata_path)
  taxonomy_df <- prep_load_taxonomy(taxonomy_path)
  genes_df <- prep_load_metagenome(genes_path)
  pathways_df <- prep_load_pathways(pathways_path)
  if(include_metabs) metabolites_df <- prep_load_metabolites(metabolites_path)
  
  # Run processing steps that are common to all datasets and independent of dataset combinations
  datasets <- list(`T` = taxonomy_df,
                   `G` = genes_df,
                   `P` = pathways_df)
  if(include_metabs) datasets$M <- metabolites_df
    
  clean_dfs <- lapply(seq_along(datasets), function(i) {
    dataset <- datasets[[i]]
    feature_set <- names(datasets)[i]

    # 1. Sort by sample ID so the folds will be the same across different feature types
    dataset <- dataset[order(dataset$`sample_id__`),]
    
    # 2. Remove rare/constant features
    dataset <- prep_sanitize_dataset(dataset, feature_set)
    
    # 3. Fix column names
    # (Some column names may be problematic for formulas, we thus fix them)
    colnames(dataset) <- make.names(colnames(dataset))
    return(dataset)
  })
  names(clean_dfs) <- names(datasets)
  
  # Run CLR transformations on relevant feature sets
  if (clr_trans) {
    ps_count <- 0.0000001
    clean_dfs$`T` <- prep_clr_transpose(clean_dfs$`T`, ps_count)
    clean_dfs$`G` <- prep_clr_transpose(clean_dfs$`G`, ps_count)
    clean_dfs$`P` <- prep_clr_transpose(clean_dfs$`P`, ps_count)
  }
  
  # Log transform metabolite values
  # (First - replace zeros/NAs with half the minimal non-zero value)
  if (include_metabs) 
    clean_dfs$M <- clean_dfs$M %>%
      mutate_if(is.numeric, function(x) ifelse(x==0, NA, x)) %>%
      mutate_if(is.numeric, function(x) ifelse(is.na(x), min(x, na.rm = T)/2, x)) %>%
      mutate_if(is.numeric, log)
  
  # Create final list of feature sets + combinations
  datasets = list(`T` = clean_dfs$`T`,
                  `G` = clean_dfs$`G`,
                  `P` = clean_dfs$`P`)
  t_and_g <- prep_join_features(clean_dfs$`T`, clean_dfs$`G`)
  t_and_p <- prep_join_features(clean_dfs$`T`, clean_dfs$`P`)
  t_g_p   <- prep_join_features(t_and_g, clean_dfs$`P`)
  
  if ("M" %in% feature_set_types) datasets$`M` <- clean_dfs$`M`
  if ("T+G" %in% feature_set_types) datasets$`T+G` <- t_and_g
  if ("T+P" %in% feature_set_types) datasets$`T+P` <- t_and_p
  if ("T+M" %in% feature_set_types) datasets$`T+M` <- prep_join_features(clean_dfs$`T`, clean_dfs$`M`)
  if ("T+G+P" %in% feature_set_types) datasets$`T+G+P` <- t_g_p
  if ("T+G+M" %in% feature_set_types) datasets$`T+G+M` <- prep_join_features(t_and_g, clean_dfs$`M`)
  if ("T+P+M" %in% feature_set_types) datasets$`T+P+M` <- prep_join_features(t_and_p, clean_dfs$`M`)
  if ("T+G+P+M" %in% feature_set_types) datasets$`T+G+P+M` <- prep_join_features(t_g_p, clean_dfs$`M`)
  
  # Process only combinations actually required for run
  datasets <- datasets[names(datasets) %in% feature_set_types]
  
  # Run these final processing steps for each data combination
  finalized_dfs <- lapply(seq_along(datasets), function(i) {
    dataset <- datasets[[i]]
    dataset_name <- names(datasets)[i]
    
    # 1. Add disease state (this will be the label we wish to predict)
    dataset <- prep_join_metadata(dataset, metadata_df)
    dataset <- dataset %>% select(-`sample_id__`)
    
    # 2. Cluster highly-redundant features 
    # (We randomly choose a single representative, and document identified clusters)
    dataset <- cluster_cluster_features(
      dataset = dataset, 
      feature_set_type = dataset_name,
      cluster_type = cluster_type,
      clusters_output = clusters_output
    )
    
    return(dataset)
  })
  names(finalized_dfs) <- names(datasets)
  
  return(finalized_dfs)
}

########################################################################
# This script includes functions for post-processing of result files
#  generated from the ml_pipeline.R script.
# The post_prepare_rdata() function saves summarized results text files 
#  and an R data file for convenience.
########################################################################

require(config)
require(metap)
require(readr)
require(dplyr)
require(stringr)

source('src/ml_pipeline/utils.R')


# Read configuration file
config <- config::get(file = "src/ml_pipeline/config.yml")


# Combine all cv_results (saved as one file per dataset) into a single table, 
#  re-format some columns and add extra information.
# Each row in the output represents a single fold/repeat in a specific dataset 
#  using a specific setting ("run_name").
# Also writes the table to a file by default.
post_load_cv_results <- function(dir_path = config$paths$results_tables_dir, 
                                 output_file = config$paths$combined_cv_results) {
  files <- list.files(dir_path, pattern = "\\_pipeline.csv$")
  all_data <- bind_rows(lapply(files, function(result_file_name) {
    df <- read_csv(file.path(dir_path, result_file_name), 
                   col_types = cols(.default = "?", feature_set_type = "c"), 
                   show_col_types = FALSE)
    return(df)
  }))
  
  # Remove the extra plus at the beginning, added by ml_pipeline.R for excel
  all_data$run_name <- substring(all_data$run_name, 2, nchar(all_data$run_name))
  all_data$run_name <- gsub("\"", "", all_data$run_name)
  all_data <- all_data %>%
    mutate(run_name = factor(run_name)) 
  
  # Write to file
  if (! is.null(output_file)) {
    utils_save_tsv_table(all_data, output_file, sep = ',')
    message('Wrote CV results to file: ', output_file)
  }
  
  return(all_data)
}


# Combine all feature importance results (saved as one file per dataset) into a single table, 
#  re-format some columns and add extra information.
# Also writes the table to a file by default.
post_load_feature_importance <- function(dir_path = config$paths$results_tables_dir, 
                                         output_file = config$paths$combined_feature_importance,
                                         cv_results = NULL) {
  
  # Read all feature importance files available in results folder
  files <- list.files(dir_path, pattern = "\\_feature_importance.csv$")
  all_data <- bind_rows(lapply(files, function(result_file_name) {
    df <- read_csv(file.path(dir_path, result_file_name), 
                   show_col_types = FALSE)
    return(df)
  }))
  
  # Remove the extra plus at the beginning, added by ml_pipeline.R for excel
  all_data$run_name <- substring(all_data$run_name, 2, nchar(all_data$run_name))
  all_data$run_name <- gsub("\"", "", all_data$run_name)
  
  # Join extra data from cv results
  if (is.null(cv_results)) cv_results <- post_load_cv_results(dir_path, output_file = NULL)
  tmp_cv_results <- cv_results %>%
    dplyr::select(dataset, run_name, mean_out_of_fold_test_auc) %>%
    distinct()
  all_data <- all_data %>%
    inner_join(tmp_cv_results, by = c('dataset', 'run_name'))
  
  # Add flag if a feature was part of a cluster
  all_data <- all_data %>% 
    mutate(is_cluster_rep = grepl("_Cluster[0-9]+$", feature)) %>%
    mutate(cluster_id = ifelse(is_cluster_rep, gsub("^.*_Cluster", "", feature), NA))
  
  # Add feature type
  all_data <- all_data %>% 
    mutate(feature_type = substr(feature, 1, 1)) %>%
    mutate(feature_set_type = gsub("^.* ", "", run_name))
  
  # Add per feature the number of times it was selected and other summary stats
  #  (min = 1, max = n_folds X n_repeats). 
  #  Features never selected are not in the list.
  feat_imp_sum <- all_data %>% 
    group_by(dataset, run_name, feature, feature_type, is_cluster_rep, cluster_id) %>% 
    summarise(n_times_selected = n(), 
              mean_importance = mean(importance),
              combined_p = post_combine_p_vals_fisher(pvalue),
              .groups = "drop") %>%
    ## Add FDR
    group_by(dataset, run_name) %>% 
    mutate(combined_fdr = p.adjust(combined_p, method = "fdr")) %>%
    ungroup()
  
  # Add these summary stats to main table
  all_data <- all_data %>%
    left_join(feat_imp_sum, 
              by = c("dataset",
                     "run_name",
                     "feature",
                     "feature_type",
                     "is_cluster_rep",
                     "cluster_id"))
  
  # Write to file
  if (! is.null(output_file)) {
    utils_save_tsv_table(all_data, output_file)
    message('Wrote feature importance results to file: ', output_file)
  }
    
  return(all_data)
}

post_summarize_feat_imp <- function(feat_imp,
                                    fdr_threshold = 0.1,
                                    n_times_selected_threshold = 0.5,
                                    n_total_folds = config$outer_n_folds * config$outer_n_repeats) {
  feat_imp_sum <- feat_imp %>%
    dplyr::select(dataset, 
                   run_name, 
                   mean_out_of_fold_test_auc, 
                   feature, 
                   feature_type,
                   feature_set_type,
                   mean_importance, 
                   combined_p, 
                   combined_fdr, 
                   n_times_selected, 
                   is_cluster_rep, 
                   cluster_id) %>%
    distinct() %>%
    mutate(contributor = 
             (combined_fdr <= fdr_threshold) & 
             (n_times_selected/n_total_folds > n_times_selected_threshold))
  return(feat_imp_sum)
}

# Prepares an rdata object with all results in it, as well as saves text files with combined results (combined over all datasets)
post_prepare_rdata <- function(dir_path = config$paths$results_tables_dir, 
                               output_file = config$paths$results_rdata, 
                               input_files_dir = config$paths$ml_input_dir,
                               auc_threshold = 0.6) {
  # 1. Load ML modelling results
  ######################################
  
  ## Load (and combine into one table) all CV results 
  cv_results <- post_load_cv_results()
  message("Completed reading all CV result files")
  
  # 2. Load feature importance results 
  ######################################
  
  feat_imp <- post_load_feature_importance(cv_results = cv_results)
  feat_imp <- feat_imp %>% filter(grepl("\\-Sh", run_name))
  message("Completed reading all feature importance files")
  feat_imp_sum <- post_summarize_feat_imp(feat_imp)
  
  # 3. Summarize statistics about each dataset
  ############################################
  
  cv_datasets_summary <- cv_results %>% 
    dplyr::select(dataset, feature_set_type,
                  tuned, shuffled,
                  fs_type, n_healthy, n_disease, 
                  n_features_origin_T, n_features_origin_S, 
                  n_features_origin_P, n_features_origin_M, 
                  mean_out_of_fold_test_auc) %>% 
    distinct() 
  
  # 4. Mark low-performing datasets
  ######################################
  
  datasets_to_analyze <- cv_results %>%
    group_by(dataset) %>%
    filter(max(mean_out_of_fold_test_auc, na.rm = TRUE) >= auc_threshold) %>%
    pull(dataset) %>%
    unique()
  
  datasets_discarded <- 
    setdiff(unique(cv_results$dataset), 
            datasets_to_analyze)
  
  message("The following datasets have AUC over ", auc_threshold, ": ",
          paste(datasets_to_analyze, collapse = ", "))
  message("Discarded datasets: ",
          paste(datasets_discarded, collapse = ", "))
  
  feat_imp <- feat_imp %>%
    filter(dataset %in% datasets_to_analyze) 
  
  feat_imp_sum <- feat_imp_sum %>%
    filter(dataset %in% datasets_to_analyze)

  # 5. Load clusters data
  ######################################
  
  all_clusters <- bind_rows(
    lapply(list.files(input_files_dir, pattern = "clusters.tsv$", recursive = T),
           function (result_file_name) {
             d <- dirname(result_file_name)
             df <- read_csv(file.path(input_files_dir, result_file_name), 
                            col_types = cols(.default = "?", feature_set_type = "c"), 
                            show_col_types = FALSE)
             df$dataset <- d
             return(df)
           })) %>%
    filter(dataset %in% datasets_to_analyze) %>%
    rename(feature_orig = feature) %>%
    rename(feature = representative)
  
  # 6. Save
  ######################################
  
  save(cv_results, 
       cv_datasets_summary,
       feat_imp, 
       feat_imp_sum, 
       all_clusters, 
       datasets_to_analyze, 
       datasets_discarded,
       file = output_file)
}


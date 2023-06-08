require(config)
require(logger)
require(tidyverse)
require(broom.helpers)
require(readr)

source("Boruta/boruta.R")

fs_boruta <- function(data_df, quantileThresh = 0.9, verbose = TRUE, maxRuns = 50, withTentative = TRUE, addLogs = TRUE) {
  logs <- list()
  
  results <- Boruta(DiseaseState ~ .,
                    data = data_df,
                    pValue = 0.01,
                    quantileThresh = quantileThresh,
                    maxRuns = maxRuns,
                    doTrace = 0)
  
  logs[['boruta_n_confirmed']] <- results$finalDecision[results$finalDecision == 'Confirmed'] %>% length
  logs[['boruta_n_rejected']] <- results$finalDecision[results$finalDecision == 'Rejected'] %>% length
  logs[['boruta_n_tentative']] <- results$finalDecision[results$finalDecision == 'Tentative'] %>% length
  
  if (verbose) 
    log_debug("Boruta: {logs[['boruta_n_confirmed']]} confirmed, {logs[['boruta_n_rejected']]} rejected, {logs[['boruta_n_tentative']]} tentative.")
  
  selected_cols <- getSelectedAttributes(results, withTentative = withTentative) %>%
    .clean_backticks()  # for variable names with spaces
  
  selected_cols <- c('DiseaseState', selected_cols)
  filtered_data <- data_df %>% select(all_of(selected_cols))

  if (! addLogs) logs <- list()
  
  return(list(
    logs = logs,
    filtered_data = filtered_data
  ))
}


fs_utest <- function(data_df,
                     thresh_type = 'pvalue',
                     thresh_pvalue = 0.25,
                     thresh_percentage = 0.2,
                     thresh_count = 20,
                     fallback = TRUE,
                     fallback_n = 20,
                     cap = TRUE,
                     cap_n = 1000) {
  logs <- list()
  
  healthy_ds <- data_df %>% filter(DiseaseState == 'healthy')
  disease_ds <- data_df %>% filter(DiseaseState == 'disease')
  
  cols <- colnames(data_df)  
  cols <- cols[cols != 'DiseaseState']
  
  p_values <- sapply(cols, function(col_name) {
    return(wilcox.test(healthy_ds[[col_name]], disease_ds[[col_name]])$p.value)
  })
  
  results <- data.frame(col=cols, p_value=unname(p_values))
  results$corrected_p_values <- p.adjust(results$p_value, method = "fdr")
  results <- results[order(results$p_value),]
  
  if (thresh_type == 'pvalue')
    filtered_results <- results %>% filter(corrected_p_values <= thresh_pvalue)
  else if (thresh_type == 'count')
    filtered_results <- results[1:thresh_count, ]
  else if (thresh_type == 'percentage')
    filtered_results <- results[1:(as.integer(nrow(results) * thresh_percentage)), ]
  else
    stop('U test: Wrong threshold type. Available options: pvalue, count, percentage')
  
  # Log
  logs[['utest_n_confirmed']] <- nrow(filtered_results)
  
  # If not enough features were selected, take the top X features (best p values) so that we can still train a model
  if (fallback & nrow(filtered_results) < fallback_n) {
    logs[['utest_fallback_or_cap']] <- 'Y'
    log_debug("Only {nrow(filtered_results)} features were selected. Fallback to {fallback_n} features.")
    filtered_results <- results[1:fallback_n, ]
  }
  
  # If an extremely high number of features were selected, take only the top X features (best p values) so that we can still train a model
  if (fallback & nrow(filtered_results) < fallback_n) {
    logs[['utest_fallback_or_cap']] <- 'Y'
    log_debug("Over {cap_n} features were selected. Taking only top {cap_n} features.")
    filtered_results <- results[1:cap_n, ]
  }
  
  selected_cols <- c('DiseaseState', filtered_results[['col']])
  filtered_data <- data_df %>% select(all_of(selected_cols))

  return(list(
    logs = logs,
    filtered_data = filtered_data
  ))
}

# The concept of ensemble feature selection: https://link.springer.com/content/pdf/10.1007/978-3-540-87481-2_21.pdf
fs_ensemble <- function(data_df, ens_method = "intersection") {
  logs <- list()
  
  # Run both u-test method and Boruta-90
  u_results <- fs_utest(data_df)
  b90_results <- fs_boruta(data_df, quantileThresh = 0.9)
  
  # Add logs from both methods
  logs <- c(logs, u_results$logs, b90_results$logs)
  
  # Now take the union/intersection of features selected
  if (ens_method == "union") {
    filtered_data <- bind_cols(
      u_results$filtered_data, 
      b90_results$filtered_data %>% select(-any_of(names(u_results$filtered_data)))
    )
    log_debug("Ensemble feature selection: {ncol(u_results$filtered_data)-1} selected by u-test, {ncol(b90_results$filtered_data)-1} selected by Boruta. Union includes {ncol(filtered_data)-1} final features.")
  } else if (ens_method == "intersection") {
    filtered_data <- u_results$filtered_data %>% 
      select(any_of(names(b90_results$filtered_data)))
    log_debug("Ensemble feature selection: {ncol(u_results$filtered_data)-1} selected by u-test, {ncol(b90_results$filtered_data)-1} selected by Boruta. Intersection includes {ncol(filtered_data)-1} final features.")
  }
  
  # Add final number of selected features to log
  logs[['ensemble_fs_n_total_selected']] <- ncol(filtered_data) - 1
  
  return(list(
    logs = logs,
    filtered_data = filtered_data
  ))
}


fs_perturb_boruta <- function(data_df, n_perturbs = 20, perc_samples = 0.8, consens_threshold = 0.9, borutaQuantileThresh = 0.9) {
  logs <- list()
  selected_feats <- list()
  
  # At each iteration choose a random subset of samples, run boruta, and record the selected features
  set.seed(777)
  for (i in 1:n_perturbs) {
    b90_results <- fs_boruta(data_df %>% slice_sample(prop = perc_samples), quantileThresh = borutaQuantileThresh, maxRuns = 50, addLogs = FALSE)
    selected_feats[[i]] <- names(b90_results$filtered_data)
  }
  
  # Now take only features selected at least <consens_threshold> of the time
  perturb_counts <- data.frame(table(unlist(selected_feats)))
  selected_cols <- c(
    'DiseaseState', 
    perturb_counts %>%
      mutate(Var1 = as.character(Var1)) %>%
      filter(Var1 != 'DiseaseState') %>%
      filter(Freq >= n_perturbs * consens_threshold) %>%
      pull(Var1)
  )
  filtered_data <- data_df %>% select(all_of(selected_cols))
  
  # Add final number of selected features to log
  logs[['boruta-perturbations_n_total_selected']] <- ncol(filtered_data) - 1
  
  return(list(
    logs = logs,
    filtered_data = filtered_data
  ))
}

# Run boruta twice or more to deal with extremely high numbers of features
fs_repeated_boruta <- function(data_df, quantileThresh = 0.9, repeats = 2) {
  logs <- list()
  selected_feats <- list()
  
  filtered_data <- data_df
  for (i in 1:repeats) {
    # Run boruta
    b90_results <- fs_boruta(filtered_data, quantileThresh = quantileThresh, addLogs = FALSE)
    # Update feature table
    filtered_data <- b90_results$filtered_data
  }
  
  # Add final number of selected features to log
  logs[['double-boruta_n_total_selected']] <- ncol(filtered_data) - 1
  
  return(list(
    logs = logs,
    filtered_data = filtered_data
  ))
}

# Select features by Altmann's feature importance
fs_altmann <- function(data_df, p_threshold = 0.1, n = NULL, addLogs = TRUE) {
  logs <- list()
  
  # Create a recipe
  data_recipe <- data_df %>%
    recipe(DiseaseState ~ .) %>%
    prep()
  
  # Create a RF object
  capped_mtry <- min(30, floor(sqrt(ncol(data_df)-1)))
  rf_obj <- 
    rand_forest(mode = "classification", mtry = capped_mtry) %>%
    set_mode("classification") %>%
    set_engine("ranger", importance = "permutation", num.threads = config::get("ranger_n_threads"))
  
  # Fit RF
  final_model <- workflow() %>%
    add_recipe(data_recipe) %>%
    add_model(rf_obj) %>% 
    fit(data = data_df)
  fit_obj <- final_model %>% extract_fit_parsnip()
  
  # Get altman's p's
  feature_importance <- ranger::importance_pvalues(fit_obj$fit,
                                                   method = 'altmann',
                                                   num.permutations = 100,
                                                   formula = as.formula(data_recipe),
                                                   data = data_df)
  
  # Take only features with p below threshold / top n features (i.e. n lowest p values)
  if (!is.null(n)) {
    feature_importance <- feature_importance[sample(nrow(feature_importance)),] # Randomize rows so that in case of p value ties no omic is favored
    feature_importance <- feature_importance[order(feature_importance[,"pvalue"]),]
    selected_feats <- rownames(feature_importance)[1:n]
  } else {
    selected_feats <- rownames(feature_importance)[feature_importance[,"pvalue"] < p_threshold]
  }
  
  selected_cols <- c(selected_feats, 'DiseaseState')
  filtered_data <- data_df %>% select(all_of(selected_cols))

  # Add final number of selected features to log
  logs[['altman_n_total_selected']] <- ncol(filtered_data) - 1
  log_debug("Altman feature selection: {logs[['altman_n_total_selected']]} features selected")
  if (! addLogs) logs <- list()
  
  return(list(
    logs = logs,
    filtered_data = filtered_data
  ))
}


# Main feature selection function
fs_select_features <- function(train_df,
                               test_df,
                               fs_type) {
  logs <- list()
  
  if (fs_type == 'none') {
    return(list(
      logs = logs,
      train_df = train_df,
      test_df = test_df
    ))
  } else if (fs_type == 'pertBoruta') {
    results <- fs_perturb_boruta(train_df)
  } else if (fs_type == 'ensemble') {
    results <- fs_ensemble(train_df)
  } else if (fs_type == 'utest') {
    results <- fs_utest(train_df)
  } else if (fs_type == 'Boruta80') {
    results <- fs_boruta(train_df, quantileThresh = 0.8)
  } else if (fs_type == 'Boruta90') {
    results <- fs_boruta(train_df, quantileThresh = 0.9)
  } else if (fs_type == 'RepeatBoruta90') {
    results <- fs_repeated_boruta(train_df, quantileThresh = 0.9)
  } else if (fs_type == 'altmann') {
    results <- fs_altmann(train_df)
  } else if (fs_type == 'altmann_top20') {
    results <- fs_altmann(train_df, n = 20)
  } else {
    stop('Wrong FS type.')
  }
  
  logs <- c(logs, results$logs)
  filtered_data <- results$filtered_data
  
  logs[['n_features_fs_selected']] <- ncol(filtered_data) - 1
  
  # If no features were selected, take all the features
  if (ncol(filtered_data) == 1) {
    filtered_data <- train_df
    log_warn("0 features selected by feature selection. Taking all features instead.")
  }
  
  # Match the test data to have the same features as the train data
  test_df <- test_df %>% select(names(filtered_data))

  return(list(
    logs = logs,
    train_df = filtered_data,
    test_df = test_df
  ))
}

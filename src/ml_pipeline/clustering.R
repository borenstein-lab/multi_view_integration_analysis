require(config)
require(logger)
require(amap)

source('src/ml_pipeline/utils.R')

cluster_calculate_clusters <- function(dataset,
                                       dist_method,
                                       max_dist) {
  # Cluster the features with the chosen distance method
  features_only <- dataset %>%
    select(-DiseaseState) %>%
    t()
  
  dist_mat <- Dist(features_only, method = dist_method)
  hclust_res <- hclust(dist_mat, method = 'complete')
  clusters <- cutree(hclust_res, h = max_dist)
  
  # Create a data frame with the features in each cluster
  clusters_df <- data.frame(feature = names(clusters), cluster = unname(clusters)) %>%
    group_by(cluster) %>%
    mutate(singleton = n() == 1) %>%
    # Randomly choose a representative per cluster
    # We mark representatives of non-singleton clusters with the cluster ID
    mutate(is_representative = feature == first(feature)) %>%
    mutate(representative = 
             ifelse(singleton,
                    first(feature),
                    paste0(first(feature),'_Cluster',cluster))) %>%
    ungroup()
  
  # Create a map from current feature names to new feature names 
  #  (a name is changed only if the feature became a representative of a non-singleton cluster)
  feature_name_map <- clusters_df$representative[clusters_df$is_representative]
  names(feature_name_map) <- clusters_df$feature[clusters_df$is_representative]
  
  return(list(clusters_df = clusters_df, feature_name_map = feature_name_map))
}

cluster_init_clusters_table <- function(clusters_output) {
  utils_save_tsv_table(data.frame("feature" = character(0),
                                  "cluster" = character(0),
                                  "singleton" = character(0),
                                  "representative" = character(0),
                                  "feature_set_type" = character(0)), 
                       clusters_output,
                       sep = ',')
}

# The function can either:
# 1. Get one dataset, cluster it and return it.
# 2. Get a train and a test datasets, cluster the train, change the test features
#    to match the clustered train features, and return both new train and test.
cluster_cluster_features <- function(dataset = NULL,
                                     train = NULL,
                                     test = NULL,
                                     feature_set_type,
                                     cluster_type,
                                     clusters_output,
                                     dist_method = 'abscorrelation') {
  logs <- list()
  is_train_test <- is.null(dataset)
  
  if (cluster_type == 'none') {
    if (!is_train_test)
      return(dataset)
    else
      return(list(train = train, test = test))
  } else if (cluster_type == 'clustering99') {
    max_dist <- 0.01
  } else if (cluster_type == 'clustering95') {
    max_dist <- 0.05
  } else {
    stop('Wrong clustering type')
  }

  dataset_for_clustering <- if (!is_train_test) dataset else train
  clusters <- cluster_calculate_clusters(dataset_for_clustering, dist_method, max_dist)
  log_debug("Clustering: {nrow(clusters$clusters_df)} features grouped into {n_distinct(clusters$clusters_df$cluster)} clusters, in dataset {feature_set_type}")
  
  # Filter only cluster representatives
  all_representatives <- clusters$clusters_df %>%
    filter(is_representative) %>%
    pull(feature) 
  relevant_features <- c(all_representatives, 'DiseaseState')
  clustered_dataset <- dataset_for_clustering %>% select(all_of(relevant_features))
  
  # Rename relevant features 
  # (features that are representatives of bigger clusters)
  feature_name_map <- clusters$feature_name_map
  feature_name_map <- c('DiseaseState' = 'DiseaseState', feature_name_map)
  clustered_dataset <- clustered_dataset %>%
    rename_with(function(x) {return(feature_name_map[x])})

  # If there is also a test dataset, apply the clustered features to it
  if (is_train_test) {
    clustered_test <- test %>%
      select(all_of(relevant_features)) %>% 
      rename_with(function(x) {return(feature_name_map[x])})
  }
  
  # Save clusters data for later analysis 
  utils_save_tsv_table(clusters$clusters_df %>%
                         select(-is_representative) %>%
                         mutate(feature_set_type = feature_set_type), 
                       clusters_output, 
                       append = TRUE,
                       col.names = FALSE,
                       sep = ',')
  
  if (!is_train_test)
    return(clustered_dataset)
  else
    return(list(train=clustered_dataset, test=clustered_test))
}

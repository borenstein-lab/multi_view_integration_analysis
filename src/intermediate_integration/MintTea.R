source("src/intermediate_integration/utils.R")

#' MintTea - A multi-view integration tool for microbiome data
#' Version 1.0
#'
#' The MintTea function receives a table containing multiple concatenated views, 
#'  alongside sample ids and study groups. It identifies modules of features 
#'  from all views that are both highly associated with one another and often 
#'  also highly associated with disease state. After identifying the modules,
#'  their association to disease is evaluated using AUROC.
#'  The function returns the modules themselves alongside overall AUROC 
#'  evaluations, module specific AUROC evaluations, and various other statistics. 
#'  The entire evaluation process is also performed on shuffled modules, 
#'  serving as null models to compare to.
#' 
#' The pipeline could be run using one or more set of pipeline parameters. If 
#'  more than one parameter is given in any of the arguments (only those 
#'  starting with the param_ prefix) then all parameter combinations will be 
#'  used, and results will include a column indicating the specific set of 
#'  parameters used for each result. User is advised to examine results and 
#'  choose a single combination of parameters that best fits his data and 
#'  research objectives.
#' 
#' Future versions will include:
#' - Support for running in parallel on multiple threads.
#' - Currently supports only "healthy" and "disease" study groups. Will be 
#'   relaxed to any arbitrary group names.
#' - Support for more than 2 study groups.
#' - Support for continuous labels.
#'
#' @param proc_data A single table containing all features of all views. Samples 
#'   are rows and features are columns. Two special columns expected are a 
#'   column holding sample identifiers and a column holding study groups 
#'   ("healthy" and "disease"). Features from each view should start with a 
#'   prefix representing the specific view followed by two underscores, e.g. 
#'   'T_' for taxonomy, 'G_' for genes, 'P_' for pathways, 'M_' for metabolites. 
#'   Features in each view should be pre-processed according to common practices. 
#'   It is advised to remove rare features.
#'   Highly correlated features should preferably be clustered (see clustering.R 
#'   script as an example).
#'   
#' @param study_group_column Name of column holding study groups 
#' 
#' @param sample_id_column Name of column holding sample identifiers
#' 
#' @param view_prefixes
#' 
#' @param param_diablo_keepX Number of features to select from each view, 
#'   serving as a constraint for the sparse CCA. Note: these are sparsity 
#'   constraints for the CCA modules, not the final consensus modules. Higher 
#'   values will produce larger modules. 
#'   More than one value can be provided if sensitivity analysis is desired.
#'   
#' @param param_diablo_design A prior on expected relations between different 
#' views. Supports values between 0 and 1 (inclusive). 0 indicates no 
#'   association between views is expected, and modules should maximize 
#'   association with disease only. 1 indicates expected inter-view association 
#'   and modules should therefore maximize both disease-association and 
#'   between-view associations.
#'   More than one value can be provided if sensitivity analysis is desired.
#'   
#' @param param_n_repeats Number of sCCA repeats on data subsamples.
#'   More than one value can be provided if sensitivity analysis is desired.
#'   
#' @param param_n_folds Number of folds used for subsetting the data before 
#'   running sCCA (DIABLO). A value of 10, for example, will divide the data 
#'   into 10 subsets, and then run CCA on 9/10 of the data, excluding each 
#'   subset one at a time. Lower values will result in smaller subsets used for 
#'   training and accordingly to higher variability between sCCA models. In 
#'   such a case we expect less modules to be identified, but their robustness 
#'   to be higher.
#'   More than one value can be provided if sensitivity analysis is desired.
#'   
#' @param param_diablo_ncomp Number of sCCA components to extract each DIABLO 
#'   run. Note that DIABLO authors recommend using only the first few components. 
#'   Typically, components >3 are less stable, and will often not contribute to 
#'   final consensus modules.
#'   More than one value can be provided if sensitivity analysis is desired.
#'   
#' @param param_edge_thresholds Number between 0 and 1 (exclusive), determining 
#'   the threshold for consensus components. Higher values mean more 
#'   conservative results. Values between 0.5 and 0.8 are recommended.
#'   More than one value can be provided if sensitivity analysis is desired.
#'   
#' @param n_evaluation_repeats Number of cross-validation repeats for overall 
#'   AUROC estimation
#'   
#' @param n_evaluation_folds Number of cross-validation folds for overall AUROC 
#'   estimation
#'   
#' @param log_level See 'library(logger); ?log_levels'
#' 
#' @param seed For result replicability
#'
#' @return A list of MintTea's multi-view modules, for each MintTea pipeline 
#'   setting used. For each module, the following properties are returned: 
#'   Module size, list of module features, module's AUROC and shuffled modules' 
#'   AUROC, average correlations between features from different views, same for 
#'   shuffled modules, and a list of edge weights in case the user wants to draw 
#'   the module as a network. 
#' 
#' If `return_main_results_only=FALSE`, then raw results are returned. 
#'   These include: 
#'   `sens_analysis_modules` - The table lists all identified modules (i.e., the full list of features in each module), for each pipeline setting.
#'   `latent_vars` - 1st prinicipal component (PC) of each module, for each pipeline setting.
#'   `module_variance_expl` - Variance explained by first PC of the true modules, as well as shuffled modules. True modules are expected to have significantly higher levels of variance explained in their first PC (compared to shuffled modules) as features are highly associated with one another.
#'   `module_inter_omic_cors` - Average correlation between features from different views, per module and per pipeline setting. Again, compared to shuffled modules.
#'   `summary_module_aucs` - AUROC of each module by itself, describing the module's association with the disease. Computed using its first PC and evaluated over repeated cross-validation. Compared to shuffled modules.
#'   `summary_overall_aucs` - Combined AUROC of all modules of a dataset, using their first PC and a simple random forest or logistic regression model, and evaluated over repeated cross-validation. Describes the overall predictive power of all modules combined. Compared to shuffled modules.
#'   
#' @export
#'
#' @examples
#'  source('src/intermediate_integration/MintTea.R')
#'  library(readr)
#'  proc_data <- read_delim("data/example_data_for_minttea/proc_data.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
#'  minttea_results <- MintTea(proc_data, view_prefixes = c('T', 'G', 'P', 'M'), param_diablo_design = c(1), param_diablo_ncomp = c(1), param_edge_thresholds = c(0.98), param_diablo_keepX = c(3))
#'  
MintTea <- function(
    proc_data, 
    study_group_column = "disease_state", 
    sample_id_column = "sample_id",
    view_prefixes,
    param_diablo_keepX = c(10),
    param_diablo_design = c(0.5),
    param_n_repeats = c(10),
    param_n_folds = c(10),
    param_diablo_ncomp = c(3),
    param_edge_thresholds = c(0.8),
    n_evaluation_repeats = 5,
    n_evaluation_folds = 10,
    log_level = 'DEBUG',
    return_main_results_only = TRUE,
    seed = 27
    ) {
  
  # Check that all required packages are installed
  installed <- rownames(installed.packages())
  if (! "dplyr" %in% installed)    stop("Please install package 'dplyr'")
  if (! "mixOmics" %in% installed) stop("Please install package 'mixOmics'")
  if (! "rsample" %in% installed)  stop("Please install package 'rsample'")
  if (! "pROC" %in% installed)     stop("Please install package 'pROC'")
  if (! "igraph" %in% installed)   stop("Please install package 'igraph'")
  if (! "ranger" %in% installed)   stop("Please install package 'ranger'")
  if (! "logger" %in% installed)   stop("Please install package 'logger'")
  if (! "conflicted" %in% installed) stop("Please install package 'conflicted'")
  if (! "stringr" %in% installed) stop("Please install package 'stringr'")
  
  # Check that 'utils.R' was sourced
  if ((!exists('is_valid_r_name')) |
      (!exists('organize_data_for_diablo')) |
      (!exists('get_auc')) ) 
    stop("Please make sure you succesfully ran `source('src/intermediate_integration/utils.R')` before executing MintTea")
  
  # Required packages
  require(dplyr)      # Tested with version: 1.0.10
  require(mixOmics)   # Tested with version: 6.18.1
  require(rsample)    # Tested with version: 1.1.0
  require(pROC)       # Tested with version: 1.18.0
  require(igraph)     # Tested with version: 1.3.4
  require(ranger)     # Tested with version: 0.14.1
  require(logger)     # Tested with version: 0.2.2
  require(conflicted) # Tested with version: 1.1.0
  require(stringr)    # Tested with version: 1.4.1
  conflict_prefer("select", "dplyr", quiet = T)
  conflict_prefer("filter", "dplyr", quiet = T)
  
  log_threshold(log_level)
  log_debug("Starting MintTea")
  
  # Parameter validations
  if (! is.numeric(param_diablo_keepX)) log_fatal('Invalid argument to *param_diablo_keepX*')
  if (min(param_diablo_keepX) < 3) log_fatal('Invalid argument to *param_diablo_keepX* (values should be >= 3)')
  if (! is.numeric(param_diablo_design)) log_fatal('Invalid argument to *param_diablo_design*')
  if (max(param_diablo_design) > 1) log_fatal('Invalid argument to *param_diablo_design* (values should be between 0 to 1)')
  if (min(param_diablo_design) < 0) log_fatal('Invalid argument to *param_diablo_design* (values should be between 0 to 1)')
  if (! is.numeric(param_n_repeats)) log_fatal('Invalid argument to *param_n_repeats*')
  if (! is.numeric(param_n_folds)) log_fatal('Invalid argument to *param_n_folds*')
  if (! is.numeric(param_diablo_ncomp)) log_fatal('Invalid argument to *param_diablo_ncomp*')
  if (min(param_diablo_ncomp) < 0) log_fatal('Invalid argument to *param_diablo_ncomp* (values should be between 1 to ~10)')
  if (max(param_diablo_ncomp) > 10) log_fatal('Invalid argument to *param_diablo_ncomp* (values should be between 1 to ~10)')
  if (! is.numeric(param_edge_thresholds)) log_fatal('Invalid argument to *param_edge_thresholds*')
  if (max(param_edge_thresholds) >= 1) log_fatal('Invalid argument to *param_edge_thresholds* (values should be between 0 to 1)')
  if (min(param_edge_thresholds) <= 0) log_fatal('Invalid argument to *param_edge_thresholds* (values should be between 0 to 1)')
  if (length(view_prefixes) > 4) log_warn('MintTea was not tested for more than 4 views.')
  if (length(view_prefixes) > 4) log_warn('MintTea was not tested for more than 4 views.')
  if (any(str_length(view_prefixes) == 0)) log_fatal('Invalid view prefix (empty character)')
  if (any(! sapply(view_prefixes, is_valid_r_name))) log_fatal('Invalid view prefix (should comply to valid R variable names. See: www.w3schools.com/r/r_variables_name.asp)')
    
  # 1. Organize input to DIABLO
  # ----------------------------------------------------------------
  diablo_input <- organize_data_for_diablo(
    proc_data, 
    study_group_column, 
    sample_id_column,
    view_prefixes,
    min_features_per_view = param_diablo_keepX)
  log_debug("Completed data preparations")
  
  # 2. Wrap DIABLO in a repeated sub-sampling procedure 
  #    to identify robust modules
  # ----------------------------------------------------------------
  
  log_debug('Running repeated sGCCA')
  
  # Function for wrapping DIABLO using specific pipeline parameters, eventually 
  #  returning a list of pairs of features and how often did they appear 
  #  together in a DIABLO component (i.e. both had non-zero loadings).
  run_wrapped_diablo <- function(
    diablo_input,
    diablo_design = 2/3,
    n_repeats = 3,
    n_folds = 10,
    diablo_keepX = 10, 
    diablo_ncomp = 5) {
    
    # Set number of variables to take from each omic
    list_keepX <- list()
    for(x in names(diablo_input$X)) list_keepX[[x]] <- rep(diablo_keepX, diablo_ncomp)
    
    # Get number of samples
    n_samples <- length(diablo_input$sample_ids)
    
    # Generate data subsets (=repeated cross-validation) 
    set.seed(seed) 
    folds <- vfold_cv(
      data.frame(sample_num = 1:n_samples), 
      v = n_folds, 
      repeats = n_repeats
    )
    
    # Place holders for AUC stats and component loadings
    cv_loadings <- data.frame()
    
    for (i in 1:nrow(folds)) {
      fold_id <- paste(folds$splits[[i]]$id, collapse = "_")
      
      # Partition the data into train and test ("held-out")
      train_samples <- folds$splits[[i]]$in_id
      held_out_samples_i <- setdiff(1:n_samples, train_samples)
      new_diablo_input <- remove_samples_from_input(diablo_input, held_out_samples_i)
      held_out_data <- remove_samples_from_input(diablo_input, train_samples)
      
      # Train DIABLO model
      tmp_diablo_out <- block.splsda(
        new_diablo_input$X, 
        new_diablo_input$Y, 
        ncomp = diablo_ncomp, 
        keepX = list_keepX, 
        max.iter = 300,
        design = diablo_design
      )
      
      # Make prediction on held out samples using DIABLO's approach 
      # predictions <- predict(tmp_diablo_out, held_out_data$X)
      #
      # Calculate auroc on held out data
      # tmp <- data.frame(
      #   preds = predictions$AveragedPredict[,1,], 
      #   true_label = held_out_data$Y
      #   )
      # cv_aucs_v <- lapply(
      #   1:diablo_ncomp, 
      #   function(i) get_auc(paste0("preds.dim", i), tmp)
      #   )
      # names(cv_aucs_v) <- paste0("auc_comp", 1:diablo_ncomp)
      # 
      # cv_aucs <- bind_rows(
      #   cv_aucs, 
      #   data.frame(cv_aucs_v) %>%
      #     mutate(fold_id = fold_id) 
      # )
      
      # Save all non-zero loadings
      tmp_loadings <- organize_diablo_loadings(tmp_diablo_out) %>%
        # Add full component identifiers (i.e. component number + specific fold)
        mutate(fold_id = fold_id) %>%
        mutate(comp_fold_id = paste0(fold_id, '_Comp', component)) %>%
        dplyr::select(-component)
      cv_loadings <- bind_rows(cv_loadings, tmp_loadings)
    }
    
    # For each pair of features count the number of times they appeared together
    # (Note - each pair currently appears twice: [A,B], [B,A])
    feat_pairs <- 
      # Get pairs of features that were found in the same component
      cv_loadings %>%
      full_join(cv_loadings, by = c("fold_id","comp_fold_id")) %>%
      filter(feature.x != feature.y) %>%
      # In cases where a pair of features come out together in more than one 
      #  component in the same run, we only count it once
      dplyr::select(feature.x, feature_set.x, feature.y, feature_set.y, fold_id) %>% 
      distinct() %>%
      group_by(feature.x, feature.y) %>%
      summarise(N = n(), .groups = "drop") 
    
    return(list(feat_pairs = feat_pairs))
  }
  
  # Run sensitivity analysis on various pipeline parameters
  feat_pairs <- data.frame()
  
  # For convenience - create a string that represents the entire parameter choices
  run_id_pattern <- 'keep_%i//des_%g//nrep_%i//nfol_%i//ncom_%i//edge_%g'
  
  for (diablo_keepX in param_diablo_keepX) {
    for (diablo_design in param_diablo_design) {
      for (n_repeats in param_n_repeats) {
        for (n_folds in param_n_folds) {
          for (diablo_ncomp in param_diablo_ncomp) {
            log_debug(paste(
              diablo_keepX,
              diablo_design,
              n_repeats,
              n_folds,
              diablo_ncomp,
              sep = ' // '
            ))
            
            # Run wrapped diablo
            tmp <- run_wrapped_diablo(
              diablo_input,
              diablo_design = diablo_design,
              n_repeats = n_repeats,
              n_folds = n_folds,
              diablo_keepX = diablo_keepX, 
              diablo_ncomp = diablo_ncomp
            )
            
            feat_pairs <- bind_rows(
              feat_pairs,
              tmp$feat_pairs %>%
                mutate(param_keepX = diablo_keepX) %>%
                mutate(param_diablo_design = diablo_design) %>%
                mutate(param_n_repeats = n_repeats) %>%
                mutate(param_n_folds = n_folds) %>%
                mutate(param_ncomp = diablo_ncomp) 
            )
          } # Done iterating over diablo_ncomp
        } # Done iterating over n_folds
      } # Done iterating over n_folds
    } # Done iterating over n_folds
  } # Done iterating over diablo_keepX
  
  # 3. Extract final components for each setting 
  #    (combination of pipeline parameters)
  # ----------------------------------------------------------------
  
  # We add one last parameter for our sensitivity analysis:
  # the 'edge_threshold' defines for a pair of 2 features, what portion of 
  # sub-sampling repeats do we require these features to be included in the same 
  # module in order for them to be placed together in our final modules.
  tmp_feat_pairs <- data.frame()
  for (p in param_edge_thresholds) {
    tmp_feat_pairs <- bind_rows(
      tmp_feat_pairs,
      feat_pairs %>% mutate(param_edge_threshold = p)
    ) 
  }
  
  # Add a "run identifier"
  tmp_feat_pairs <- tmp_feat_pairs %>%
    mutate(run_id = sprintf(
      run_id_pattern, 
      param_keepX, 
      param_diablo_design, 
      param_n_repeats, 
      param_n_folds, 
      param_ncomp, 
      param_edge_threshold
    ))
  # table(tmp_feat_pairs$run_id)
  
  # Thin table to record all runs for sensitivity analysis
  sens_analysis_runs <- tmp_feat_pairs %>%
    select(run_id, starts_with('param_')) %>%
    distinct() %>%
    tidyr::pivot_longer(
      cols = starts_with('param_'), 
      names_to = 'param', 
      names_prefix = 'param_', 
      values_to = 'value'
    )
  
  # For each setting, infer modules and organize them in a list
  log_debug('Identifying consensus multi-view modules')
  sens_analysis_modules <- data.frame()
  for (curr_run in unique(sens_analysis_runs$run_id)) {
    # Get number of overall repeats
    n_rep_folds <- sens_analysis_runs %>%
      filter(run_id == curr_run) %>%
      filter(param %in% c('n_folds','n_repeats')) %>%
      pull(value)
    n_rep_folds <- n_rep_folds[1] * n_rep_folds[2]
    
    # Get edge threshold
    edge_threshold <- sens_analysis_runs %>%
      filter(run_id == curr_run) %>%
      filter(param == 'edge_threshold') %>%
      pull(value)
    
    # Create list of edges by taking all feature pairs that appeared together 
    #  over <edge_threshold> of the times
    edges <- tmp_feat_pairs %>% 
      filter(run_id == curr_run) %>%
      filter(N >= n_rep_folds * edge_threshold)
    
    if(nrow(edges) == 0) {
      cat('\n')
      log_debug('No modules identified for MintTea setting: {curr_run}')
      next
    }
    
    vertices <- edges %>% 
      select(feature.x) %>% 
      distinct() 
    
    net <- graph_from_data_frame(
      d = edges, 
      vertices = vertices, 
      directed = F
    ) 
    
    tmp <- components(net, mode = "strong")
    final_modules <- data.frame(
      feature = names(tmp$membership),
      module = paste0('module',unname(tmp$membership)),
      run_id = curr_run
    ) %>% tibble::remove_rownames()
    
    sens_analysis_modules <- bind_rows(
      sens_analysis_modules,
      final_modules
    )
  }
  
  # Discard settings for which no modules were identified
  sens_analysis_runs <- sens_analysis_runs %>%
    filter(run_id %in% unique(sens_analysis_modules$run_id))
  
  if (nrow(sens_analysis_runs) == 0) return(list())
  
  # 4. Compute module summaries (1st PC), per setting
  # ----------------------------------------------------------------
  
  log_debug('Computing new variables based on modules, and shuffled versions for comparison')
  latent_vars <- list()
  module_variance_expl <- list()
  module_inter_omic_cors <- list()
  shuffled_modules <- list()
  set.seed(seed)
  
  # For each setting, get modules first PCs
  for (curr_run in unique(sens_analysis_runs$run_id)) {
    
    # Initialize
    latent_vars[[curr_run]] <- list()
    module_variance_expl[[curr_run]] <- list()
    module_inter_omic_cors[[curr_run]] <- list()
    shuffled_modules[[curr_run]] <- list()
    
    # # Organize current settings 
    # tmp <- strsplit(strsplit(curr_run, "//")[[1]], '_')
    # selected_settings <- lapply(tmp, function(x) as.numeric(x[2]))
    # names(selected_settings) <- sapply(tmp, function(x) x[1])
    
    # Get modules list
    selected_modules <- sens_analysis_modules %>%
      filter(run_id == curr_run) %>%
      select(-run_id)
    log_debug('With settings {curr_run}, identified {n_distinct(selected_modules$module)} modules, with {selected_modules %>% group_by(module) %>% summarise(N = n()) %>% pull(N) %>% mean() %>% round(2)} features on average')
    
    # 0 = true modules, 1-99 = shuffled versions
    # log_debug('Fetching shuffled modules for null distributions')
    for (i in 0:99) {
      shuf <- ifelse(i == 0, FALSE, TRUE)
      run <- ifelse(shuf, paste0('shuf', i), 'true')
      
      latent_vars[[curr_run]][[run]] <- data.frame(label = diablo_input$Y) # Place holder
      module_variance_expl[[curr_run]][[run]] <- data.frame()
      module_inter_omic_cors[[curr_run]][[run]] <- data.frame()
      shuffled_modules[[curr_run]][[run]] <- list()
      
      # Iterate over modules 
      for (module_id in 1:n_distinct(selected_modules$module)) {
        
        # Get module
        tmp <- selected_modules %>%
          filter(module == paste0('module', module_id))
        
        # Extract relevant data to use as PCA input
        # (take only features in this module from each omic)
        xlist <- lapply(diablo_input$X, function(x) {x[,colnames(x) %in% tmp$feature,drop=FALSE]})
        
        # Remove views without features in this module
        rem <- sapply(xlist, ncol) == 0
        xlist <- xlist[!rem]
        
        # Shuffle if needed 
        if (shuf) {
          # For each view, sample same amount of features as in true module
          xlist <- lapply(
            names(diablo_input$X),
            function(fs) {
              num_feats_to_sample <- sum(grepl(paste0("^",fs,"__"), tmp$feature))
              x <- diablo_input$X[[fs]]
              x[,sample(colnames(x), num_feats_to_sample),drop=FALSE]
            }
          )
          names(xlist) <- names(diablo_input$X)
          shuffled_modules[[curr_run]][[run]][[module_id]] <- unname(unlist(lapply(xlist, colnames)))
        }
        
        # Run PCA on features in module
        module_raw_data <- bind_cols(lapply(xlist, as.data.frame))
        pca <- prcomp(module_raw_data, center = TRUE, scale. = TRUE) 
        var_explained <- pca$sdev^2/sum(pca$sdev^2) * 100 # Percent of variance explained
        
        # Record variance explained (we expect high values for 1st components)
        module_variance_expl[[curr_run]][[run]] <- bind_rows(
          module_variance_expl[[curr_run]][[run]],
          data.frame(module = module_id,
                     var_explained_pc1 = var_explained[1],
                     var_explained_pc2 = var_explained[2]) %>% 
            tibble::remove_rownames()
        ) 
        
        # Record mean correlations between features from different views
        tmp_inter_omic_cors <- get_all_inter_omic_corrs(module_raw_data)
        module_inter_omic_cors[[curr_run]][[run]] <- bind_rows(
          module_inter_omic_cors[[curr_run]][[run]],
          data.frame(module = module_id,
                     n_pairs = nrow(tmp_inter_omic_cors),
                     avg_pears_corr = mean(tmp_inter_omic_cors$pearson_corr),
                     avg_spear_corr = mean(tmp_inter_omic_cors$spearman_corr))
        ) 
        
        # Extract and save new latent feature
        new_latent <- data.frame(pca$x[,'PC1'])
        names(new_latent) <- paste0("module", module_id)
        latent_vars[[curr_run]][[run]] <- bind_cols(
          latent_vars[[curr_run]][[run]], 
          new_latent
        )
      }
    }
    
    # Re-organize
    module_variance_expl[[curr_run]] <- bind_rows(module_variance_expl[[curr_run]], .id = 'run')
    module_inter_omic_cors[[curr_run]] <- bind_rows(module_inter_omic_cors[[curr_run]], .id = 'run')
  }
  
  
  # 5. Run repeated cross validation to assess AUC of single modules 
  #    (1st PCs) and overall (All 1st PCs of all modules)
  # ----------------------------------------------------------------
  
  # First, get the number of repeats/folds from ML results
  log_debug('Starting model evaluation using cross-validation')
  n_samples <- length(diablo_input$sample_ids)
  
  # Now generate folds consistent with ML pipeline
  set.seed(seed)
  folds <- vfold_cv(
    data.frame(sample_num = 1:n_samples),
    v = n_evaluation_folds,
    repeats = n_evaluation_repeats
  )
  
  # Collect AUC stats here
  cv_overall_auc <- data.frame()
  module_aucs <- data.frame()
  
  # Iterate over settings
  for (curr_run in unique(sens_analysis_runs$run_id)) {
    log_debug('Training and evaluating modules using settings: {curr_run}')
    
    # Get actual module list
    selected_modules <- sens_analysis_modules %>%
      filter(run_id == curr_run) %>%
      select(-run_id)
    
    # Iterate over true and shuffled modules
    for (run in names(shuffled_modules[[curr_run]])) {
      if (log_level == 'DEBUG') cat('.')
      
      # Cross validation loop
      for (i in 1:nrow(folds)) {
        fold_id <- paste(folds$splits[[i]]$id, collapse = " ")
        
        # Get train samples
        train_samples <- folds$splits[[i]]$in_id
        train_data <- data.frame(label = diablo_input$Y[train_samples])
        test_data <- data.frame(label = diablo_input$Y[-train_samples])
        
        # Iterate over modules
        for (module_id in 1:n_distinct(selected_modules$module)) {
          # Get only features in module
          tmp <- selected_modules %>%
            filter(module == paste0('module', module_id)) %>%
            pull(feature)
          
          # If we're running on shuffled modules - take shuffled modules features (override real features)
          if (run != 'true') {
            tmp <- shuffled_modules[[curr_run]][[run]][[module_id]]
          }
          module_raw_data <- lapply(diablo_input$X, function(x) {x[,colnames(x) %in% tmp, drop=FALSE]})
          
          # Remove omics without features in this module
          rem <- sapply(module_raw_data, ncol) == 0
          module_raw_data <- bind_cols(lapply(module_raw_data[!rem], as.data.frame))
          
          # Train PCA
          pca <- prcomp(module_raw_data[train_samples,], center = TRUE, scale. = TRUE) 
          
          # Run same PCA loadings on test data (test data gets scaled according to train data)
          predicted_pcs <- predict(pca, newdata = module_raw_data[-train_samples,])
          
          # Organize train + test tables
          tmp_train_data <- pca$x[,'PC1'] %>% data.frame(); names(tmp_train_data)[1] <- paste0('module',module_id)
          train_data <- bind_cols(train_data, tmp_train_data)
          tmp_test_data <- predicted_pcs[,'PC1'] %>% data.frame(); names(tmp_test_data)[1] <- paste0('module',module_id)
          test_data <- bind_cols(test_data, tmp_test_data)
          
          # Calculate and save module-specific AUC
          module_aucs <- bind_rows(
            module_aucs,
            data.frame(
              setting = curr_run,
              run = run,
              fold_id = fold_id,
              module_id = module_id,
              auc = get_auc(paste0('module', module_id), data = test_data, label_col = "label", auto_direction = TRUE)
            )
          )
        } # Done iterating over modules
        
        # Important for GLM
        train_data <- train_data %>%
          mutate(label = factor(label, levels = c("healthy", "disease")))
        test_data <- test_data %>%
          mutate(label = factor(label, levels = c("healthy", "disease")))
        
        # Train logistic regression / RF
        logistic_model <- glm(label ~ ., data = train_data, family = "binomial")
        rf <- ranger(label ~ ., data = train_data, probability = TRUE)
        
        # Make prediction on held out samples
        preds_glm <- predict(logistic_model, test_data, type = "response")
        preds_rf <- predict(rf, test_data)$predictions
        
        # Calculate and save overall auroc
        tmp <- data.frame(preds_glm = preds_glm, preds_rf = preds_rf[,'disease'], true_label = test_data$label)
        cv_overall_auc <- bind_rows(
          cv_overall_auc,
          data.frame(
            setting = curr_run,
            run = run,
            fold_id = fold_id,
            test_auc_glm = get_auc('preds_glm', tmp),
            test_auc_rf = get_auc('preds_rf', tmp)
          )
        )
      } # Completed cross-validation
    } # Completed iterations over true and shuffled modules
    if (log_level == 'DEBUG') cat('\n')
  } # Completed iterations over pipeline settings
  
  # Summarize
  summary_overall_aucs <- cv_overall_auc %>%
    group_by(setting, run) %>%
    summarise(mean_overall_rf_auc = mean(test_auc_rf, na.rm = TRUE),
              sd_overall_rf_auc = sd(test_auc_rf, na.rm = TRUE),
              mean_overall_lm_auc = mean(test_auc_glm, na.rm = TRUE),
              sd_overall_lm_auc = sd(test_auc_glm, na.rm = TRUE),
              .groups = "drop")
  
  summary_module_aucs <- module_aucs %>%
    group_by(setting, run, module_id) %>%
    summarise(mean_module_auc = mean(auc, na.rm = TRUE),
              sd_module_auc = sd(auc, na.rm = TRUE),
              .groups = "drop")
  log_debug("Completed MintTea")
  
  # 6. Return a list of identified modules or summary tables
  # ----------------------------------------------------------------
  
  # Re-organize results in tables into a list of MintTea modules
  modules_list <- list()
  # Iterate over MintTea settings (for cases when MintTea was executed with 
  #  several optional settings)
  for (curr_run in unique(sens_analysis_modules$run_id)) {
    modules_list[[curr_run]] <- list()
    
    # Iterate over modules identified 
    for (module_id in (sens_analysis_modules %>% 
                       filter(run_id == curr_run) %>% 
                       pull(module) %>% 
                       unique())) {
      
      # Extract list of features of this module
      module_feats <- sens_analysis_modules %>% 
        filter(run_id == curr_run) %>% 
        filter(module == module_id) %>%
        pull(feature)
      module_id_num <- gsub('module','',module_id) %>% as.integer()
        
      # Populate module properties and evaluation statistics 
      #  (extract info from the various result tables)
      modules_list[[curr_run]][[module_id]] <- list()
      modules_list[[curr_run]][[module_id]][['module_size']] <- length(module_feats)
      modules_list[[curr_run]][[module_id]][['features']] <- module_feats
      modules_list[[curr_run]][[module_id]][['module_edges']] <- 
        tmp_feat_pairs %>% 
          filter(run_id == curr_run) %>%
          filter((feature.x %in% module_feats) & (feature.y %in% module_feats)) %>%
          filter(feature.x < feature.y) %>%
          select(feature.x, feature.y, N) %>%
          mutate(edge_weight = round(N / (n_repeats * n_folds),4))
      modules_list[[curr_run]][[module_id]][['auroc']] <- 
        summary_module_aucs %>%
          filter(setting == curr_run) %>% 
          filter(module_id == module_id_num) %>% 
          filter(run == 'true') %>% 
          pull(mean_module_auc)
      modules_list[[curr_run]][[module_id]][['shuffled_auroc']] <- 
        summary_module_aucs %>%
          filter(setting == curr_run) %>% 
          filter(module_id == module_id_num) %>% 
          filter(run != 'true') %>% 
          pull(mean_module_auc)
      modules_list[[curr_run]][[module_id]][['inter_view_corr']] <- 
        module_inter_omic_cors[[curr_run]] %>% 
          filter(module == module_id_num) %>% 
          filter(run == 'true') %>% 
          pull(avg_pears_corr)
      modules_list[[curr_run]][[module_id]][['shuffled_inter_view_corr']] <- 
        module_inter_omic_cors[[curr_run]] %>% 
          filter(module == module_id_num) %>% 
          filter(run != 'true') %>% 
          pull(avg_pears_corr)
    }
  }
  
  if (return_main_results_only)
    return(modules_list)
    
  return(list(diablo_input = diablo_input, 
              feat_pairs = feat_pairs,
              sens_analysis_modules = sens_analysis_modules,
              sens_analysis_runs = sens_analysis_runs,
              latent_vars = latent_vars,
              module_variance_expl = module_variance_expl,
              module_inter_omic_cors = module_inter_omic_cors,
              module_aucs = module_aucs,
              cv_overall_auc = cv_overall_auc,
              summary_module_aucs = summary_module_aucs,
              summary_overall_aucs = summary_overall_aucs))
}

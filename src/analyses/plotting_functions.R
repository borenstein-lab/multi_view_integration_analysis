###################################################################
# Plotting functions for manuscript figures and other explorations
###################################################################

config <- config::get(file = "src/ml_pipeline/config.yml")

# Gets a list of p values and returns significance stars
get_signif_mark <- function(ps) {
  sapply(ps, function(p) {
    if (p < 0.0001) return("***")
    if (p < 0.001) return("**")
    if (p < 0.05) return("*")
    if (p < 0.1) return(".")
    return("")
  })
}

# Get informative names for features
get_pwy_names <- function(pwy_codes) {
  pwy_codes <- make.names(pwy_codes)
  
  # MetaCyc pathways
  pwy_names <- read_delim(
    config$paths$metacyc_pathways2taxa, 
    delim = "\t", 
    escape_double = FALSE, 
    trim_ws = TRUE, 
    col_select = c(Pathways, `Common-Name`),
    show_col_types = FALSE) %>%
    rename(ID = 1, NAME = 2)
  
  # Kegg modules
  kegg_mod_names <- read_delim(
    config$paths$kegg_module_names, 
    delim = "\t", 
    escape_double = FALSE, 
    trim_ws = TRUE, 
    show_col_types = FALSE) %>%
    rename(ID = 1, NAME = 2)
  pwy_names <- bind_rows(pwy_names, kegg_mod_names)
  pwy_names$ID_FIXED <- make.names(pwy_names$ID)
  
  pwy_map <- pwy_names$NAME
  names(pwy_map) <- pwy_names$ID_FIXED
  
  pwys_with_no_name <- pwy_codes[! pwy_codes %in% names(pwy_map)]
  if (length(pwys_with_no_name) > 0)
    stop(paste('The following pathway codes could not be mapped to names:',
               paste(pwys_with_no_name, collapse = ', ')))
  
  return(unname(pwy_map[pwy_codes]))
}


get_feature_descriptions <- function(features) {
  features_uniq <- unique(features)
  
  T_features <- grep('^T__', features_uniq, value = TRUE)
  P_features <- grep('^P__', features_uniq, value = TRUE)
  M_features <- grep('^M__', features_uniq, value = TRUE)
  S_features <- grep('^S__', features_uniq, value = TRUE)
  
  T_features_map <- gsub('\\.', ' ', gsub('(^T__|_Cluster[0-9]+$)', '', T_features))
  P_features_map <- gsub('(^P__|_Cluster[0-9]+$)', '', P_features) %>% get_pwy_names()
  M_features_map <- gsub('(^M__|_Cluster[0-9]+$)', '', M_features)
  S_features_map <- gsub('(^S__|_Cluster[0-9]+$)', '', S_features) 
  
  P_features_map <- gsub('(<sub>|</sub>|<i>|</i>|<I>|</I>)', '_', P_features_map)
  
  T_clustered <- grepl('_Cluster[0-9]+$', T_features)
  P_clustered <- grepl('_Cluster[0-9]+$', P_features)
  M_clustered <- grepl('_Cluster[0-9]+$', M_features)
  S_clustered <- grepl('_Cluster[0-9]+$', S_features)
  
  T_features_map <- paste0(T_features_map, ifelse(T_clustered, ' [Clustered features]', ''))
  P_features_map <- paste0(P_features_map, ifelse(P_clustered, ' [Clustered features]', ''))
  M_features_map <- paste0(M_features_map, ifelse(M_clustered, ' [Clustered features]', ''))
  S_features_map <- paste0(S_features_map, ifelse(S_clustered, ' [Clustered features]', ''))
  
  names(T_features_map) <- T_features
  names(P_features_map) <- P_features
  names(M_features_map) <- M_features
  names(S_features_map) <- S_features
  
  feature_name_map <- c(T_features_map, P_features_map, M_features_map, S_features_map)
  
  return(unname(feature_name_map[features]))
}

get_u_tests <- function(df, diablo_input) {
  all_utests <- data.frame()
  for (d in unique(df$dataset)) {
    feats_to_test <- df$feature[df$dataset == d]
    proc_data <- bind_cols(lapply(diablo_input[[d]]$X, as.data.frame))
    dis_labels <- diablo_input[[d]]$Y
    for (f in feats_to_test) {
      h_vec <- proc_data[[f]][dis_labels == "healthy"]
      d_vec <- proc_data[[f]][dis_labels != "healthy"]
      utest <- wilcox.test(h_vec, d_vec, digits.rank = 7, exact = FALSE)
      
      all_utests <- bind_rows(
        all_utests,
        data.frame(
          dataset = d,
          feature = f,
          p = utest$p.value,
          mean_abun_healthy = mean(h_vec), 
          mean_abun_disease = mean(d_vec)
        )
      )
    }
  }
  
  all_utests <- all_utests %>%
    group_by(dataset) %>%
    mutate(fdr_utest = p.adjust(p, method = "fdr")) %>%
    ungroup() %>%
    mutate(increased_in = ifelse(mean_abun_healthy > mean_abun_disease, "healthy", "disease")) %>%
    mutate(increased_in = ifelse(fdr_utest > 0.05, "non-significant", increased_in)) %>%
    select(-p)
  return(all_utests)
}

# Plot basic statistics about datasets and ml pipeline results
plot_basic_stats <- function(cv_results,  
                             dataset_order = NULL, 
                             feature_type_color_map) {
  
  # Extract needed columns
  tmp <- cv_results %>%
    select(dataset, n_healthy, n_disease, 
           run_name, shuffled, feature_set_type, fold_id,  
           n_features_origin_T, n_features_origin_S, 
           n_features_origin_P, n_features_origin_M, 
           n_features_for_train_final_T, n_features_for_train_final_S, 
           n_features_for_train_final_P, n_features_for_train_final_M, 
           out_of_fold_test_auc, mean_out_of_fold_test_auc) 
  
  # Reorder datasets
  if (!is.null(dataset_order)) {
    tmp <- tmp %>% 
      filter(dataset %in% dataset_order) %>%
      mutate(dataset = factor(dataset, levels = dataset_order[dataset_order %in% tmp$dataset]))
  } else { 
    tmp <- tmp %>% 
      mutate(dataset = factor(dataset))
  }
  
  # ---------------- Strip 1: sample size per dataset ----------------->
  tmp2 <- tmp %>%
    filter(! shuffled) %>%
    group_by(dataset, n_healthy, n_disease) %>%
    summarise(n_feature_sets = n_distinct(feature_set_type), .groups = 'drop') %>%
    mutate(label = paste0(" C: ", n_healthy,
                          "\n D: ", n_disease)) %>%
    arrange(dataset)
  
  p1_col_width = 0.8
  
  p1 <- ggplot(tmp2, aes(x = dataset, y = n_healthy + n_disease)) +
    geom_rect(data = tmp2 %>% select(dataset) %>% distinct(), 
              aes(alpha = dataset %>% as.numeric() %% 2 == 0), 
              xmin = -Inf, xmax = Inf, ymin = -Inf, 
              ymax = Inf, fill = 'gray80', inherit.aes = FALSE) +
    scale_alpha_manual(values = c('FALSE' = 0, 'TRUE' = 0.3), guide = "none") +
    geom_col(color="black", fill = "grey80", width = p1_col_width) +
    geom_text(aes(label=label), vjust = 0.5, hjust = 0, size = 2.9, lineheight = .9) +
    geom_col(aes(y = n_healthy), color = "black", fill = "grey50", width = p1_col_width) +
    scale_y_continuous(expand = c(0, 0, 0.48, 0)) +
    coord_flip() +
    theme_classic() +
    xlab(NULL) +
    ylab("Number of\nsamples") +
    facet_grid(rows = vars(dataset), space = "free_y", scales = "free_y", switch = "y") +
    theme(plot.margin = unit(c(5.5,6.5,7.8,5.5), "points")) +
    theme(axis.text.y = element_text(size = 10)) +
    theme(panel.grid.major.x = element_line(linewidth = 0.5, color = "grey93")) +
    theme(axis.title.x = element_text(size = 11)) +
    theme(strip.background = element_blank()) +
    theme(strip.placement = "outside") +
    theme(strip.text.y.left = element_text(size = 10, angle = 0, hjust = 1)) +
    theme(panel.spacing.y = unit(6, "points")) +
    theme(axis.text.y.left = element_blank())
  
  # Patch to make figure heights same as following plots
  # gt <- ggplot_gtable(ggplot_build(p1))
  # patch_coefs <- (gt$heights[7] * nrow(tmp2)) * (tmp2$n_feature_sets / sum(tmp2$n_feature_sets))
  # for (i in 1:nrow(tmp2)) {
  #   gt$heights[5 + 2*i] <- patch_coefs[i]
  # }
  # p1 <- ggplotify::as.ggplot(gt)
  
  # ---- Strip 2: N features per dataset per view ---->
  tmp2 <- tmp %>%
    filter(! shuffled) %>%
    select(dataset, feature_set_type,
           starts_with('n_features_orig')) %>%
    distinct() %>%
    tidyr::pivot_longer(
      cols = starts_with('n_features_orig'), 
      names_to = 'feature_type', 
      names_prefix = 'n_features_origin_', 
      values_to = 'n_features'
    ) %>%
    mutate(feature_type = factor(feature_type, levels = rev(c('T','P','M','S')))) %>%
    filter(n_features > 0)
  
  p2 <- ggplot(tmp2, 
               aes(x = feature_type, 
                   y = n_features, 
                   fill = feature_type)) +
    geom_col(color="black", alpha = 0.2, width = 0.8, position = "stack") + # Plotted twice as patch
    geom_rect(data = tmp2 %>% select(dataset) %>% distinct(), 
              aes(alpha = dataset %>% as.numeric() %% 2 == 0), 
              xmin = -Inf, xmax = Inf, ymin = -Inf, 
              ymax = Inf, fill = 'gray80', inherit.aes = FALSE) +
    scale_alpha_manual(values = c('FALSE' = 0, 'TRUE' = 0.3), guide = "none") +
    geom_col(color="black", alpha = 0.8, width = 0.8, position = "stack") +
    scale_x_discrete(expand = c(0, 1.1)) +
    scale_y_continuous(expand = c(0, 0, 0.1, 0)) +
    coord_flip() +
    theme_classic() +
    facet_grid(rows = vars(dataset), scales = "free_y", switch = "y") +
    xlab(NULL) +
    ylab("Number of\nfeatures") +
    scale_fill_manual(name = "View", values = feature_type_color_map,
                      guide = guide_legend(reverse = TRUE)) +
    theme(panel.grid.major.y = element_blank()) +
    theme(panel.grid.major.x = element_line(linewidth = 0.5, color = "grey93")) +
    theme(axis.text.y = element_blank()) +
    # theme(legend.position = "none") +
    theme(axis.title.x = element_text(size = 11)) +
    theme(strip.background = element_blank(), strip.text = element_blank()) +
    theme(panel.spacing.y = unit(6, "points"))
  
  tmp_rel_widths = c(3.7,3)
  plot_grid(p1, p2, 
            nrow = 1, 
            rel_widths = tmp_rel_widths, 
            align = 'h', axis = 'tb')
}

# Plot basic statistics about datasets and ml pipeline results
plot_basic_stats2 <- function(cv_results, feature_type_color_map) {
  
  tmp <- cv_results %>%
    filter(! shuffled) %>%
    mutate(n_samples = n_healthy+n_disease) %>%
    select(dataset, 
           n_samples,
           feature_set_type,
           starts_with('n_features_orig')) %>%
    select(-n_features_origin) %>%
    distinct() %>%
    mutate(dataset = factor(dataset, labels = paste0(dataset,'\nN = ',n_samples))) %>%
    tidyr::pivot_longer(
      cols = starts_with('n_features_orig'), 
      names_to = 'feature_type', 
      names_prefix = 'n_features_origin_', 
      values_to = 'n_features'
    ) %>%
    mutate(feature_type = factor(feature_type, levels = rev(c('T','P','M','S')))) %>%
    filter(n_features > 0)
  
  p <- ggplot(tmp, 
               aes(x = feature_type, 
                   y = n_features, 
                   fill = feature_type)) +
    geom_col(color="black", alpha = 0.2, width = 0.8, position = "stack") + # Plotted twice as patch
    geom_rect(data = tmp %>% select(dataset) %>% distinct(), 
              aes(alpha = dataset %>% as.numeric() %% 2 == 0), 
              xmin = -Inf, 
              xmax = Inf, 
              ymin = -Inf, 
              ymax = Inf, 
              fill = 'gray80', 
              inherit.aes = FALSE) +
    scale_alpha_manual(values = c('FALSE' = 0, 'TRUE' = 0.3), guide = "none") +
    geom_col(color="black", alpha = 0.8, width = 0.8, position = "stack") +
    scale_x_discrete(expand = c(0, 1.1)) +
    scale_y_continuous(expand = c(0, 0, 0.05, 0)) +
    coord_flip() +
    theme_classic() +
    facet_grid(rows = vars(dataset), scales = "free_y", switch = "y") +
    xlab(NULL) +
    ylab("Number of\nfeatures") +
    scale_fill_manual(name = NULL, values = feature_type_color_map,
                      guide = guide_legend(reverse = TRUE)) +
    theme(panel.grid.major.y = element_blank()) +
    theme(panel.grid.major.x = element_line(linewidth = 0.5, color = "grey93")) +
    theme(axis.text.y = element_blank()) +
    # theme(legend.position = "none") +
    theme(axis.title.x = element_text(size = 11)) +
    #theme(strip.background = element_blank(), strip.text = element_blank()) +
    theme(strip.background = element_blank()) +
    theme(strip.placement = "outside") +
    theme(strip.text.y.left = element_text(size = 10, angle = 0, hjust = 1)) +
    theme(panel.spacing.y = unit(6, "points")) +
    theme(legend.position = 'bottom')
  
  print(p)
}

# Same style as basic_stats plot, this time describing the diablo components in each dataset
plot_module_stats <- function(sens_analysis_modules, 
                              modules_overview, 
                              feature_type_color_map,
                              dataset_order,
                              show_rf = TRUE,
                              hide_y_axis_text = FALSE) {
  
  # ---- Strip 1: N features per dataset per module ---->
  tmp1 <- sens_analysis_modules %>%
    mutate(feature_type = substr(feature,1,1)) %>%
    mutate(module = as.numeric(gsub('module', '', module))) %>%
    group_by(dataset, module, feature_type) %>%
    summarise(N = n(), .groups = "drop") 
  
  # Only plot modules with at least 2 types of features 
  modules_to_plot <- modules_overview %>% 
    filter(multi_view) %>%
    select(dataset, module, is_interesting) %>%
    mutate(module = as.numeric(gsub('module','',module)))
 
  tmp1 <- inner_join(tmp1, modules_to_plot, by = c('dataset','module'))
     
  tmp1$module2 <- factor(paste('Module', tmp1$module),
                            levels = paste('Module', max(tmp1$module):1))
  
  p1 <- ggplot(tmp1 %>%
                 filter(dataset %in% dataset_order) %>%
                 mutate(dataset = factor(dataset, levels = dataset_order)), 
               aes(x = module2, y = N, fill = feature_type)) +
    geom_rect(data = tmp1 %>% 
                filter(dataset %in% dataset_order) %>% 
                mutate(dataset = factor(dataset, levels = dataset_order)) %>%
                select(dataset) %>% 
                distinct(),
              aes(alpha = dataset %>% as.numeric() %% 2 == 0),
              xmin = -Inf, 
              xmax = Inf, 
              ymin = -Inf,
              ymax = Inf, 
              fill = 'gray80', 
              inherit.aes = FALSE) +
    scale_alpha_manual(values = c('FALSE' = 0, 'TRUE' = 0.3), guide = "none") +
    new_scale("alpha") +
    geom_bar(aes(alpha = is_interesting, color = is_interesting),
             width = 0.8, 
             position="stack", 
             stat="identity") +
    scale_alpha_manual(values = c('FALSE' = 0.3, 'TRUE' = 0.9), guide = "none") +
    scale_color_manual(values = c('FALSE' = 'gray70', 'TRUE' = 'black'), guide = "none") +
    geom_hline(yintercept = 0) +
    scale_y_continuous(expand = c(0, 0, 0.1, 0)) +
    scale_x_discrete(expand = c(0, 0.5)) +
    coord_flip() +
    theme_classic() +
    facet_grid(rows = vars(dataset), space = "free_y", scales = "free_y", switch = "y") +
    xlab(NULL) +
    ylab("No. of features\nof each type") +
    scale_fill_manual(name = "View", values = feature_type_color_map) +
    theme(panel.grid.major.x = 
            element_line(linewidth = 0.5, color = "grey93")) +
    theme(panel.grid.minor.x = element_blank()) +
    theme(panel.grid.major.y = 
            element_line(linewidth = 0.5, color = "grey93")) +
    theme(axis.text.y = element_text(size = 9)) +
    theme(legend.position = "none") +
    theme(axis.title.x = element_text(size = 11)) +
    theme(strip.background = element_blank()) +
    theme(strip.placement = "outside") +
    theme(strip.text.y.left = element_text(size = 10, angle = 0, hjust = 1)) +
    theme(panel.spacing.y = unit(6, "points")) 
  
  if (hide_y_axis_text) p1 <- p1 + theme(axis.text.y = element_blank())
    
  # ---- Strip 3: AUC ---->
  tmp2 <- modules_overview %>% 
    select(-is_interesting) %>%
    mutate(module = as.numeric(gsub('module','',module))) %>%
    inner_join(modules_to_plot, by = c('dataset','module'))
  
  tmp2$module2 <- factor(paste('Module', tmp2$module),
                            levels = levels(tmp1$module2))
  points_size <- ifelse(hide_y_axis_text, 3, 3.7)
  
  p2 <- ggplot(tmp2 %>%
                 filter(dataset %in% dataset_order) %>%
                 mutate(dataset = factor(dataset, levels = dataset_order)), 
               aes(x = module2)) +
    geom_rect(data = tmp1 %>% 
                filter(dataset %in% dataset_order) %>% 
                mutate(dataset = factor(dataset, levels = dataset_order)) %>%
                select(dataset) %>% 
                distinct(), 
              aes(alpha = dataset %>% as.numeric() %% 2 == 0), 
              xmin = -Inf, 
              xmax = Inf, 
              ymin = -Inf, 
              ymax = Inf, 
              fill = 'gray80', 
              inherit.aes = FALSE) +
    scale_alpha_manual(values = c('FALSE' = 0, 'TRUE' = 0.3), guide = "none") +
    geom_hline(yintercept = 0.5, color = "darkred", linetype = "dashed", linewidth = 1) +
    # Shuffled standard deviations
    geom_linerange(aes(ymax = mean_module_auc_shuffled + sd_module_auc_shuffled, 
                       ymin = mean_module_auc_shuffled - sd_module_auc_shuffled), 
                   alpha = 0.4, linewidth = 2, color = "grey70") +
    # Shuffled AUCs
    geom_point(aes(y = mean_module_auc_shuffled), 
               shape = 16, size = points_size - 0.5, 
               color = 'grey60', alpha = 0.8) +
    # True standard deviations
    # geom_linerange(aes(ymax = sdev_high_, ymin = sdev_low_, color = feature_set_type),
    #                position = position_dodge(width = p3_dodging), 
    #                alpha = 0.35, linewidth = 2) + 
    # True AUCs
    geom_point(aes(y = mean_module_auc, fill = is_interesting, color = is_interesting),
               shape = 23, 
               size = points_size, alpha = 0.9) +
    scale_fill_manual(values = c('FALSE' = '#D7E7ED', 'TRUE' = 'skyblue4'), guide = "none") +
    scale_color_manual(values = c('FALSE' = 'gray60', 'TRUE' = 'black'), guide = "none") +
    scale_y_continuous(breaks = seq(0.5,1,0.1), expand = expansion(mult = c(0.05,0.1))) +
    scale_x_discrete(expand = c(0, 0.5)) +
    coord_flip() +
    theme_classic() +
    xlab(NULL) +
    ylab("Module AUC") +
    facet_grid(rows = vars(dataset), space = "free_y", scales = "free_y", switch = "y") +
    theme(panel.grid.major.x = 
            element_line(linewidth = 0.5, color = "grey93")) +
    theme(panel.grid.major.y = 
            element_line(linewidth = 0.5, color = "grey93")) +
    theme(axis.title.x = element_text(size = 11)) +
    theme(axis.text.y = element_blank()) +
    theme(strip.background = element_blank(), strip.text = element_blank()) +
    theme(panel.spacing.y = unit(6, "points"))
  
  if (show_rf) p2 <- p2 + geom_hline(aes(yintercept = mean_auc_rf), color = "goldenrod2", linewidth = 2, alpha = 0.7)
  
  # ---- Strip 2: Cross-view correlations ---->
  p3 <- ggplot(tmp2 %>%
                 filter(dataset %in% dataset_order) %>%
                 mutate(dataset = factor(dataset, levels = dataset_order)), 
               aes(x = module2)) + 
    geom_rect(data = tmp1 %>% 
                filter(dataset %in% dataset_order) %>% 
                mutate(dataset = factor(dataset, levels = dataset_order)) %>%
                select(dataset) %>% 
                distinct(), 
              aes(alpha = dataset %>% as.numeric() %% 2 == 0), 
              xmin = -Inf, xmax = Inf, ymin = -Inf, 
              ymax = Inf, fill = 'gray80', inherit.aes = FALSE) +
    scale_alpha_manual(values = c('FALSE' = 0, 'TRUE' = 0.3), guide = "none") +
    # Shuffled standard deviations
    geom_linerange(aes(ymax = avg_spear_corr_shuffled + sd_spear_corr_shuffled, 
                       ymin = avg_spear_corr_shuffled - sd_spear_corr_shuffled), 
                   alpha = 0.4, linewidth = 2, color = "grey70") +
    # Shuffled correlations
    geom_point(aes(y = avg_spear_corr_shuffled), 
               shape = 16, size = points_size - 0.5, color = 'grey60', alpha = 0.8) +
    geom_point(aes(y = avg_spear_corr, fill = is_interesting, color = is_interesting),
               shape = 23, 
               size = points_size, 
               alpha = 0.9) +
    scale_fill_manual(values = c('FALSE' = '#D9B9AB', 'TRUE' = 'sienna4'), guide = "none") +
    scale_color_manual(values = c('FALSE' = 'gray60', 'TRUE' = 'black'), guide = "none") +
    scale_x_discrete(expand = c(0, 0.5)) +
    scale_y_continuous(expand = expansion(mult = c(0.05,0.1))) +
    coord_flip() +
    theme_classic() +
    xlab(NULL) +
    ylab("Cross-view avg.\ncorrelation") +
    facet_grid(rows = vars(dataset), space = "free_y", scales = "free_y", switch = "y") +
    theme(panel.grid.major.x = 
            element_line(linewidth = 0.5, color = "grey93")) +
    theme(panel.grid.major.y = 
            element_line(linewidth = 0.5, color = "grey93")) +
    theme(axis.title.x = element_text(size = 11)) +
    theme(axis.text.y = element_blank()) +
    theme(strip.background = element_blank(), strip.text = element_blank()) +
    theme(panel.spacing.y = unit(6, "points"))
  
  # Combine plots
  tmp_rel_widths <- c(8.1, 3, 3)
  if (hide_y_axis_text) tmp_rel_widths <- c(7, 3, 3)
  print(plot_grid(p1, p3, p2, 
                  nrow = 1, 
                  rel_widths = tmp_rel_widths, 
                  align = 'h', axis = 'tb'))
}

plot_overall_modules_aucs <- function(rf_results, summary_aucs, datasets_to_focus_on) {
  
  caption_early_integration <- "Early integration\n(Average no. of features after feature selection)"
  caption_minttea <- "Intermediate integration with MintTea\n(No. of features [No. of modules])"
  caption_minttea_shuf <- "Intermediate integration with MintTea - Shuffled"
  
  tmp <- bind_rows(
    rf_results %>%
      mutate(label_n_features = as.character(round(mean_n_features_for_train, 1))) %>%
      select(dataset, mean_auc, sd_auc, label_n_features) %>%
      mutate(model = 'rf') %>%
      mutate(pipeline = caption_early_integration) %>%
      mutate(shuf = FALSE),
    summary_aucs %>%
      filter(run == 'true') %>%
      rename(mean_auc = mean_overall_rf_auc, sd_auc = sd_overall_rf_auc) %>%
      mutate(label_n_features = paste0(n_features, ' [', n_modules, ']')) %>%
      select(dataset, mean_auc, sd_auc, label_n_features) %>%
      mutate(model = 'rf') %>%
      mutate(pipeline = caption_minttea) %>%
      mutate(shuf = FALSE),
    summary_aucs %>%
      filter(run != 'true') %>%
      group_by(dataset) %>%
      summarize(mean_auc = mean(mean_overall_rf_auc), sd_auc = sqrt(wtd.var(mean_overall_rf_auc^2))) %>%
      mutate(model = 'rf') %>%
      mutate(pipeline = caption_minttea_shuf) %>%
      mutate(shuf = TRUE)
  ) %>%
    filter(dataset %in% datasets_to_focus_on) %>%
    mutate(dataset = factor(dataset, levels = rev(datasets_to_focus_on))) %>%
    mutate(error_low = mean_auc-sd_auc) %>%
    mutate(error_high = pmin(1, mean_auc+sd_auc)) %>%
    mutate(pipeline = factor(pipeline, 
                             levels = c(caption_minttea_shuf,
                                        caption_minttea,
                                        caption_early_integration)))
  
  p_dodging = 0.8
  
  coloring <- c('lightgrey','goldenrod', 'darkred')
  names(coloring) <- c(caption_minttea_shuf,caption_minttea,caption_early_integration)
  
  p <- ggplot(tmp %>% filter(model == 'rf'), 
              mapping = aes(x = dataset, group = pipeline)) +
    # Long grey lines
    geom_linerange(aes(ymax = mean_auc), 
                   ymin = 0, color = 'darkgrey',
                   position = position_dodge(width = p_dodging), 
                   alpha = 0.6, linewidth = 0.6) + 
    geom_hline(yintercept = 0.5, color = "darkred", linetype = "dashed", linewidth = 1) +
    # Standard errors
    geom_linerange(aes(ymax = error_high, ymin = error_low),
                   position = position_dodge(width = p_dodging), 
                   alpha = 0.4, linewidth = 2, color = "grey70") +
    # AUCs
    geom_point(aes(y = mean_auc, fill = pipeline), 
               color = "black", size = 3, shape = 21,
               alpha = 0.85, position = position_dodge(width = p_dodging)) +
    # Number of features
    geom_text(aes(label = label_n_features, y = mean_auc+0.05), 
              hjust = 0, size = 3, vjust = 0.3,
              position = position_dodge(width = p_dodging)) +
    scale_y_continuous(breaks = seq(0.5,1,0.1), expand = expansion(mult = c(0, 0.15), add = 0)) +
    scale_x_discrete(expand = c(0, 0.5)) +
    scale_fill_manual(name = "Pipeline", values = coloring) +
    guides(fill = guide_legend(reverse=TRUE)) +
    coord_flip() +
    theme_classic() +
    xlab(NULL) +
    ylab("Random forest AUC") +
    theme(panel.grid.major.x = 
            element_line(linewidth = 0.5, color = "grey93")) +
    theme(legend.spacing.y = unit(1.0, 'cm'))  +
    theme(axis.title.x = element_text(size = 11)) 
  print(p)
}

plot_module_sensitivity_analysis <- function(sens_analysis_modules, latent_vars, base_module, selected_settings, settings_order = NULL) {
  require(pROC)
  
  # Final set of features for a given module
  base_features <- sens_analysis_modules %>% 
    filter(dataset == base_module$d) %>%
    filter(module == paste0('module', base_module$mod_id)) %>%
    filter(run_id == selected_settings) %>%
    pull(feature)
  
  # All relevant components/features
  df <- sens_analysis_modules %>% 
    filter(dataset == base_module$d) %>%
    group_by(run_id, module) %>%
    filter(sum(feature %in% base_features) > 0) %>%
    ungroup() 
  
  if (is.null(settings_order)) {
    # Make the selected setting first in order
    tmp_settings <- unique(sens_analysis_modules$run_id)
    settings_order <- c(selected_settings, tmp_settings[tmp_settings != selected_settings])
  }
  df <- df %>% 
    mutate(run_id = factor(
      run_id, 
      levels = settings_order
      ))
  
  # Hack for nice ordering of features
  feats_order <- sort(base_features) 
  for (i in 2:length(settings_order)) {
    tmp_feats <- df %>% 
      filter(run_id == settings_order[i]) %>%
      pull(feature)
    feats_order <- c(feats_order, tmp_feats[! tmp_feats %in% feats_order])
  }
  df <- df %>%
    mutate(feature = factor(feature, levels = rev(feats_order)))
  
  # Control module colors
  module_color <- df %>%
    select(run_id, module) %>%
    distinct() %>%
    group_by(run_id) %>%
    mutate(module_color = as.character(row_number())) %>%
    ungroup()
  
  # Add AUC's per module
  module_color$AUC <- NA
  
  # For each module plotted...
  for (i in 1:nrow(module_color)) {
    # Get PC1's of these modules
    mod_id <- module_color$module[i]
    curr_run_id <- as.character(module_color$run_id[i])
    tmp <- latent_vars[[base_module$d]][[curr_run_id]]$true
    # Get ROC's
    tmp_roc <- roc(tmp$label, tmp[[mod_id]], levels = c('healthy', 'disease'), quiet = TRUE)
    # Save
    module_color$AUC[i] <- round(tmp_roc$auc, 3)
  }
  
  # Add new properties to data to plot
  df <- df %>%
    left_join(module_color, by = c("module", "run_id"))
  
  settings_order_labels <- gsub('// nfol', '//\nnfol', gsub('//', ' // ', settings_order))
  
  df <- df %>% 
    mutate(run_id = factor(
      run_id, 
      levels = settings_order, 
      labels = settings_order_labels
    ))
  
  p1 <- ggplot(df, aes(y = feature, x = run_id, fill = module_color)) +
    geom_point(size = 2.5, shape = 21, color = 'black') +
    annotate("rect", xmin = 0.6, xmax = 1.6, ymin = 0, ymax = 1+n_distinct(df$feature), alpha = .1, fill = "cadetblue4") +
    scale_x_discrete(position = 'top') +
    scale_fill_manual(values = c('1' = '#9AD2CB', 
                                 '2' = '#E0D68A', 
                                 '3' = '#8E443D', 
                                 '4' = '#511730')) +
    theme_classic() +
    xlab(NULL) +
    ylab(NULL) +
    theme(legend.position = 'none') +
    theme(axis.text.x.top = element_text(angle = 60, hjust = 0)) +
    theme(axis.text.y = element_text(size = 6)) +
    theme(plot.margin = unit(c(10,95,1,1), 'points')) 
  
  p2 <- ggplot(module_color %>%
                 tidyr::complete(run_id, module_color) %>%
                 mutate(run_id = factor(
                   run_id, 
                   levels = settings_order, 
                   labels = settings_order_labels
                 )) %>%
                 mutate(AUC = ifelse(is.na(AUC), 0.501, AUC)), 
               aes(y = AUC, x = run_id, fill = module_color, group = module_color)) +
    geom_col(width = 0.7, color = 'black', position = 'dodge') +
    annotate("rect", xmin = 0.5, xmax = 1.5, ymin = 0, ymax = 1, alpha = .1, fill = "cadetblue4") +
    coord_cartesian(ylim = c(0.5, 1), expand = T) +
    scale_fill_manual(values = c('1' = '#9AD2CB', '2' = '#E0D68A', '3' = '#8E443D', '4' = '#511730')) +
    scale_x_discrete(expand = expansion(add = 0.7)) +
    theme_classic() +
    ylab('Module\'s\n1st PC AUC') +
    xlab(NULL) +
    theme(legend.position = 'none') +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    theme(axis.title.y = element_text(size = 10)) +
    theme(plot.margin = unit(c(10,95,1,1), 'points')) 
  
  return(list(p1 = p1, p2 = p2))
}


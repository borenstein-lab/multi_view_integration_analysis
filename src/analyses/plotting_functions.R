###################################################################
# Plotting functions for manuscript figures and other explorations
###################################################################

# Adjust dataset names for final plots if needed
fix_dataset_name <- function(df, datasets_names_to_fix) {
  for (i in 1:length(datasets_names_to_fix)) {
    df <- df %>%
      mutate(dataset = ifelse(dataset == names(datasets_names_to_fix)[i], 
                              unname(datasets_names_to_fix[i]), 
                              dataset))
  }
  return(df)
}

# Get informative names for features
get_pwy_names <- function(pwy_codes) {
  pwy_codes <- make.names(pwy_codes)
  
  pwy_names <- read_delim(
    config::get('paths')$metacyc_pathways2taxa, 
    delim = "\t", 
    escape_double = FALSE, 
    trim_ws = TRUE, 
    col_select = c(Pathways, `Common-Name`),
    show_col_types = FALSE)
  
  pwy_names$id_fixed <- make.names(pwy_names$Pathways)
  pwy_map <- pwy_names$`Common-Name`
  names(pwy_map) <- pwy_names$id_fixed
  
  pwys_with_no_name <- pwy_codes[! pwy_codes %in% names(pwy_map)]
  if (length(pwys_with_no_name) > 0)
    stop(paste('The following pathway codes could not be mapped to names:',
               paste(pwys_with_no_name, collapse = ', ')))
  
  return(unname(pwy_map[pwy_codes]))
}

get_ko_names <- function(ko_codes) {
  ko_names <- read.table(
    file = config::get('paths')$ko_names,
    header = TRUE,
    quote="", 
    fill = TRUE,
    sep = "\t")
  
  ko_map <- ko_names$DEFINITION
  names(ko_map) <- make.names(ko_names$ENTRY)
  
  # Patch
  ko_map <- c(ko_map, c('None' = 'None'))
  
  kos_with_no_name <- ko_codes[! ko_codes %in% names(ko_map)]
  if (length(kos_with_no_name) > 0)
    stop(paste('The following KO codes could not be mapped to names:',
               paste(kos_with_no_name, collapse = ', ')))
  
  return(unname(ko_map[ko_codes]))
}

get_feature_descriptions <- function(features) {
  features_uniq <- unique(features)
  
  T_features <- grep('^T__', features_uniq, value = TRUE)
  P_features <- grep('^P__', features_uniq, value = TRUE)
  M_features <- grep('^M__', features_uniq, value = TRUE)
  G_features <- grep('^G__', features_uniq, value = TRUE)
  
  T_features_map <- gsub('\\.', ' ', gsub('(^T__|_Cluster[0-9]+$)', '', T_features))
  P_features_map <- gsub('(^P__|_Cluster[0-9]+$)', '', P_features) %>% get_pwy_names()
  M_features_map <- gsub('(^M__|_Cluster[0-9]+$)', '', M_features)
  G_features_map <- gsub('(^G__|_Cluster[0-9]+$)', '', G_features) %>% get_ko_names()
  
  P_features_map <- gsub('(<sub>|</sub>|<i>|</i>|<I>|</I>)', '_', P_features_map)
  
  T_clustered <- grepl('_Cluster[0-9]+$', T_features)
  P_clustered <- grepl('_Cluster[0-9]+$', P_features)
  M_clustered <- grepl('_Cluster[0-9]+$', M_features)
  G_clustered <- grepl('_Cluster[0-9]+$', G_features)
  
  T_features_map <- paste0(T_features_map, ifelse(T_clustered, ' [Clustered features]', ''))
  P_features_map <- paste0(P_features_map, ifelse(P_clustered, ' [Clustered features]', ''))
  M_features_map <- paste0(M_features_map, ifelse(M_clustered, ' [Clustered features]', ''))
  G_features_map <- paste0(G_features_map, ifelse(G_clustered, ' [Clustered features]', ''))
  
  names(T_features_map) <- T_features
  names(P_features_map) <- P_features
  names(M_features_map) <- M_features
  names(G_features_map) <- G_features
  
  feature_name_map <- c(T_features_map, P_features_map, M_features_map, G_features_map)
  
  return(unname(feature_name_map[features]))
}

# Plot basic statistics about datasets and ml pipeline results
plot_basic_stats <- function(cv_results, mtg_type = "Shotgun", dataset_order = NULL, feature_type_color_map, p3_width = 4) {
  
  # Extract needed columns
  tmp <- cv_results %>%
    filter(metagenomics_type == mtg_type) %>%
    select(dataset, data_source, n_healthy, n_disease, 
           run_name, shuffled, feature_set_type, fold_id,  
           n_features_origin_T, n_features_origin_G, 
           n_features_origin_P, n_features_origin_M, 
           n_features_for_train_final_T, n_features_for_train_final_G, 
           n_features_for_train_final_P, n_features_for_train_final_M, 
           out_of_fold_test_auc, mean_out_of_fold_test_auc) %>%
    mutate(feature_set_type = forcats::fct_collapse(feature_set_type, 'T+G+P / T+G+P+M' = c('T+G+P', 'T+G+P+M')))
  
  # Update feature-set colors mapping
  feature_type_color_map['T+G+P / T+G+P+M'] = 'grey30'
  
  # Reorder datasets
  if (!is.null(dataset_order))
    tmp <- tmp %>% 
      filter(dataset %in% dataset_order) %>%
      mutate(dataset = factor(dataset, levels = dataset_order[dataset_order %in% tmp$dataset]))
  
  # ---------------- Strip 1: sample size per dataset ----------------->
  tmp2 <- tmp %>%
    filter(! shuffled) %>%
    group_by(dataset, data_source, n_healthy, n_disease) %>%
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
    scale_y_continuous(expand = c(0, 0, 0.45, 0)) +
    coord_flip() +
    theme_classic() +
    xlab(NULL) +
    ylab("No. samples") +
    facet_grid(rows = vars(dataset), space = "free_y", scales = "free_y", switch = "y") +
    theme(plot.margin = unit(c(5.5,6.5,7.8,5.5), "points")) +
    theme(axis.text.y = element_text(size = 10)) +
    theme(panel.grid.major.x = 
            element_line(linewidth = 0.5, color = "grey93")) +
    theme(axis.title.x = element_text(size = 11)) +
    theme(strip.background = element_blank()) +
    theme(strip.placement = "outside") +
    theme(strip.text.y.left = element_text(size = 10, angle = 0, hjust = 1)) +
    theme(panel.spacing.y = unit(6, "points")) +
    theme(axis.text.y.left = element_blank())
  
  # Patch to make figure heights same as following plots
  gt <- ggplot_gtable(ggplot_build(p1))
  patch_coefs <- (gt$heights[7] * nrow(tmp2)) * (tmp2$n_feature_sets / sum(tmp2$n_feature_sets))
  for (i in 1:nrow(tmp2)) {
    gt$heights[5 + 2*i] <- patch_coefs[i]
  }
  p1 <- ggplotify::as.ggplot(gt)
  
  # ---- Strip 2: N features per dataset per view ---->
  tmp2 <- tmp %>%
    filter(! shuffled) %>%
    select(dataset, data_source, feature_set_type,
           n_features_for_train_final_T, n_features_for_train_final_G, 
           n_features_for_train_final_P, n_features_for_train_final_M) %>%
    tidyr::pivot_longer(
      cols = c('n_features_for_train_final_T', 
               'n_features_for_train_final_G', 
               'n_features_for_train_final_P', 
               'n_features_for_train_final_M'), 
      names_to = 'feature_type', 
      names_prefix = 'n_features_for_train_final_', 
      values_to = 'n_features'
    ) %>%
    mutate(feature_type = factor(feature_type, levels = rev(c('T','G','P','M')))) %>%
    group_by(dataset, data_source, feature_set_type, feature_type) %>%
    summarise(average_n_features = mean(n_features, na.rm = TRUE), N = n(),
              .groups = "drop") %>%
    filter(average_n_features > 0) %>%
    mutate(feature_set_type = factor(feature_set_type, levels = rev(levels(feature_set_type)))) %>% 
    group_by(dataset, feature_set_type) %>% 
    summarise(average_n_features = sum(average_n_features), .groups = 'drop')
    
    p2 <- ggplot(tmp2, 
                 aes(x = feature_set_type, 
                     y = average_n_features, 
                     fill = feature_set_type)) +
      geom_col(color="black", alpha = 0.2, width = 0.8, position = "stack") + # Plotted twice as patch
      geom_rect(data = tmp2 %>% select(dataset) %>% distinct(), 
                aes(alpha = dataset %>% as.numeric() %% 2 == 0), 
                xmin = -Inf, xmax = Inf, ymin = -Inf, 
                ymax = Inf, fill = 'gray80', inherit.aes = FALSE) +
      scale_alpha_manual(values = c('FALSE' = 0, 'TRUE' = 0.3), guide = "none") +
      geom_col(color="black", alpha = 0.5, width = 0.8, position = "stack") +
      scale_x_discrete(expand = c(0, 1.1)) +
      scale_y_continuous(trans = log10_trans(),
                         breaks = trans_breaks("log10", function(x) 10^x),
                         labels = trans_format("log10", math_format(10^.x)),
                         expand = c(0, 0, 0.1, 0)) +
      coord_flip() +
      theme_classic() +
      facet_grid(rows = vars(dataset), space = "free_y", scales = "free_y", switch = "y") +
      xlab(NULL) +
      ylab("No. features after FS") +
      scale_fill_manual(name = "View", values = feature_type_color_map) +
      theme(panel.grid.major.x = 
              element_line(linewidth = 0.5, color = "grey93")) +
      theme(panel.grid.major.y = 
              element_line(linewidth = 0.5, color = "grey93")) +
      theme(panel.grid.minor.x = element_blank()) +
      theme(axis.text.y = element_blank()) +
      theme(legend.position = "none") +
      theme(axis.title.x = element_text(size = 11)) +
      theme(axis.text.x = element_text(size = 8)) +
      theme(strip.background = element_blank(), strip.text = element_blank()) +
      theme(panel.spacing.y = unit(6, "points"))
  
  # Patch (to allow alignment with p1)
  p2 = ggplotify::as.ggplot(ggplot_gtable(ggplot_build(p2)))
  
  # ---- Strip 3: Single-view AUC + SD ---->
  tmp2 <- tmp %>%
    group_by(dataset, shuffled, feature_set_type, mean_out_of_fold_test_auc) %>%
    summarise(n = n(), 
              sdev = sd(out_of_fold_test_auc, na.rm = TRUE), 
              .groups = "drop") %>%
    mutate(sdev_low = mean_out_of_fold_test_auc - sdev) %>%
    mutate(sdev_high = pmin(1, mean_out_of_fold_test_auc + sdev)) %>%
    mutate(feature_set_type = factor(feature_set_type, levels = rev(levels(feature_set_type)))) %>% 
    mutate(shuffled = ifelse(shuffled, 'shuf', '')) %>%
    tidyr::pivot_wider(id_cols = c(dataset, feature_set_type), 
                       names_from = shuffled, 
                       values_from = c(mean_out_of_fold_test_auc,n,sdev_low,sdev_high))
  
  p3_min_x <- 0.4
  
  p3 <- ggplot(data = tmp2, 
               mapping = aes(x = feature_set_type)) +
    geom_rect(data = tmp2 %>% select(dataset) %>% distinct(), 
              aes(alpha = dataset %>% as.numeric() %% 2 == 0), 
              xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              fill = 'gray80', inherit.aes = FALSE) +
    scale_alpha_manual(values = c('FALSE' = 0, 'TRUE' = 0.3), guide = "none") +
    geom_hline(yintercept = 0.5, color = "darkred", linetype = "dashed", linewidth = 1) +
    # Shuffled standard deviations
    geom_linerange(aes(ymax = sdev_high_shuf, 
                       ymin = pmax(p3_min_x, sdev_low_shuf)), 
                   alpha = 0.35, linewidth = 2, color = "grey70") +
    # Shuffled AUCs
    geom_point(aes(y = mean_out_of_fold_test_auc_shuf), 
               shape = 16, size = 2, color = 'grey60', alpha = 0.7) +
    # True standard deviations
    geom_linerange(aes(ymax = sdev_high_, 
                       ymin = pmax(p3_min_x, sdev_low_), 
                       color = feature_set_type), 
                   alpha = 0.35, linewidth = 2) + 
    # True AUCs
    geom_point(aes(fill = feature_set_type,
                   y = mean_out_of_fold_test_auc_), 
               shape = 23, color = "black", size = 3.7, alpha = 0.85) +
    scale_y_continuous(breaks = seq(0.5,1,0.1), limits = c(p3_min_x, NA), expand = c(0, 0, 0.04, 0)) +
    scale_x_discrete(expand = c(0, 1.1)) +
    coord_flip() +
    scale_fill_manual(name = "View", values = feature_type_color_map[rev(levels(tmp2$feature_set_type))], 
                      breaks = names(feature_type_color_map)[names(feature_type_color_map) %in% tmp2$feature_set_type]) +
    scale_color_manual(name = "View", values = feature_type_color_map[rev(levels(tmp2$feature_set_type))], 
                       breaks = names(feature_type_color_map)[names(feature_type_color_map) %in% tmp2$feature_set_type]) +
    facet_grid(rows = vars(dataset), space = "free_y", scales = "free_y", switch = "y") +
    theme_classic() +
    xlab(NULL) +
    ylab("Random forest AUC") +
    theme(panel.grid.major.x = 
            element_line(linewidth = 0.5, color = "grey93")) +
    theme(panel.grid.major.y = 
            element_line(linewidth = 0.5, color = "grey93")) +
    theme(axis.text.y = element_blank()) +
    theme(axis.title.x = element_text(size = 11)) +
    theme(strip.background = element_blank(), strip.text = element_blank()) +
    theme(panel.spacing.y = unit(6, "points"))
  
  # Patch (to allow alignment with p1)
  p3 = ggplotify::as.ggplot(ggplot_gtable(ggplot_build(p3)))
  
  tmp_rel_widths = c(3.5,2.5,p3_width)
  plot_grid(p1, p2, p3, 
            nrow = 1, 
            rel_widths = tmp_rel_widths, 
            align = 'h', axis = 'tb')
}

# Plots an overview of feature importance given by different ML pipelines
plot_feat_importance_overview <- function(
    d, 
    full_feature_importance, 
    # Formatting options
    w_p1 = 8, w_legends = 4, 
    # As this graph can be quite long, we limit the number of features 
    #  we include in the visualization using the following cutoffs:
    # Top features from each view that we want to include in visualization
    rank_cutoff_single_omic = 5, 
    # Top features from multi-view models we want to include
    rank_cutoff_multi_omic = 15,
    specific_features = NULL,
    limits_fill_p2_p3_p4 = rep(NA, 6)
) {
  
  # Select dataset
  tmp <- full_feature_importance %>%
    filter(dataset == d)
  
  # Select features to plot 
  if (is.null(specific_features)) {
    tmp <- tmp %>%
      # Top X from each pipeline
      filter(single_omic_rank <= rank_cutoff_single_omic | multi_omic_rank <= rank_cutoff_multi_omic) 
  } else {
    tmp <- tmp %>%
      # User specified
      filter(pretty_name %in% specific_features) 
  }
  
  # Mark cluster representatives
  tmp <- tmp %>% 
    mutate(
      feature_orig2 = 
        paste0(
          pretty_name,
          ifelse((!is.na(feature_rep_single_omic) & 
                    grepl('_Cluster[0-9]+$', feature_rep_single_omic)) |
                   (!is.na(feature_rep_multi_omic) & 
                      grepl('_Cluster[0-9]+$', feature_rep_multi_omic)),
                 '*', '')))
  
  # Organize data for plotting
  tmp <- tmp %>%
    mutate(fdr_utest_for_plot = -log10(fdr_utest)) %>%
    mutate(single_omic_fdr_for_plot = -log10(single_omic_fdr)) %>%
    mutate(multi_omic_fdr_for_plot = -log10(multi_omic_fdr)) %>%
    mutate(feature_type = factor(feature_type, levels = c('T','G','P','M'))) %>%
    mutate(feature_orig2 = factor(feature_orig2, levels = tmp %>% arrange(desc(single_omic_fdr)) %>% pull(feature_orig2)))
  
  p1 <- ggplot(tmp, aes(x = '', y = feature_orig2, fill = increased_in)) +
    geom_point(color = 'black', shape = 21, size = 4) +
    ylab(NULL) +
    xlab(NULL) +
    scale_x_discrete(position = "top") +
    scale_fill_manual(values = c("healthy" = "dodgerblue4", "disease" = "firebrick1", "non-significant" = "lightgrey"),
                      name = 'Increased in...') +
    theme_minimal() +
    facet_grid(feature_type ~ ., scales = "free_y", space = "free") +
    theme(strip.background = element_blank(), 
          strip.text.x = element_blank(), 
          strip.text = element_blank()) +
    theme(legend.text = element_text(size=11)) +
    theme(plot.margin = margin(5.5,5,5.5,0)) +
    theme(panel.grid = element_blank()) +
    theme(legend.key.width= unit(8, 'points'))
  
  max_tick_legend <- max(3, round(max(tmp$single_omic_fdr_for_plot, na.rm = TRUE)))
  min_tick_legend <- min(3, round(min(tmp$fdr_utest_for_plot, na.rm = TRUE)))
  p3 <- ggplot(tmp, aes(x = '', y = feature_orig2, 
                        fill = single_omic_fdr_for_plot)) +
    geom_point(color = 'black', shape = 22, size = 7.5) +
    geom_text(aes(label = single_omic_rank), size = 2.8) +
    ylab(NULL) +
    xlab('Single-view\nfeature importance') +
    scale_x_discrete(position = "top") +
    scale_fill_gradient(high = 'darkolivegreen4', 
                        low = lighten('darkolivegreen1', amount = 0.8),
                        na.value = "grey85",
                        name = 'Altmann FDR\nSingle-view model',
                        #breaks = seq(from = min_tick_legend, to = max_tick_legend, by = round((max_tick_legend-min_tick_legend)/4)),
                        labels = function(x) parse(text = paste0(" ~ 10^-", x)),
                        limits = limits_fill_p2_p3_p4[3:4],
                        guide = guide_colorbar(frame.colour = "black", 
                                               ticks.colour = "black")) +
    theme_minimal() +
    facet_grid(feature_type~., scales = "free_y", space = "free") +
    theme(strip.background = element_blank(), 
          strip.text.x = element_blank(), 
          strip.text = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(axis.title.x.top = element_text(size = 10, angle = 90, hjust = 0, vjust = 0.5, margin = margin(0,0,0,0))) +
    theme(legend.text = element_text(size=10)) +
    theme(plot.margin = margin(5.5,0,5.5,0)) +
    theme(panel.grid = element_blank()) +
    theme(legend.key.width= unit(8, 'points'), legend.key.height = unit(15, 'points'))
  
  max_tick_legend <- max(3, round(max(tmp$multi_omic_fdr_for_plot, na.rm = TRUE)))
  min_tick_legend <- min(3, round(min(tmp$fdr_utest_for_plot, na.rm = TRUE)))
  p4 <- ggplot(tmp %>% mutate(highlight = (multi_omic_rank <= rank_cutoff_multi_omic)), 
               aes(x = '', y = feature_orig2, fill = multi_omic_fdr_for_plot)) +
    geom_point(color = 'black', shape = 22, size = 7.5) +
    geom_text(aes(label = multi_omic_rank, color = highlight, fontface = highlight), size = 2.8) +
    scale_color_manual(values = c('FALSE' = 'grey50', 'TRUE' = 'black')) +
    scale_discrete_manual('fontface', values = c('FALSE' = 'plain', 'TRUE' = 'bold')) +
    guides(color = 'none', fontface = 'none') +
    ylab(NULL) +
    xlab('Multi-view\nfeature importance') +
    scale_x_discrete(position = "top") +
    scale_fill_gradient(high = 'darkorange3', low = 'blanchedalmond',
                        na.value = "grey85",
                        name = 'Altmann FDR\nMulti-view model',
                        #breaks = seq(from = min_tick_legend, to = max_tick_legend, by = round((max_tick_legend-min_tick_legend)/4)),
                        labels = function(x) parse(text = paste0(" ~ 10^-", x)),
                        limits = limits_fill_p2_p3_p4[5:6],
                        guide = guide_colorbar(frame.colour = "black", 
                                               ticks.colour = "black")) +
    theme_minimal() +
    facet_grid(feature_type~., scales = "free_y", space = "free") +
    theme(strip.background = element_blank(), 
          strip.text.x = element_blank(), 
          strip.text = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(axis.title.x.top = element_text(size = 10, angle = 90, hjust = 0, vjust = 0.5, margin = margin(0,0,0,0))) +
    theme(plot.margin = margin(5.5,0,5.5,0)) +
    theme(legend.text = element_text(size = 10)) +
    theme(panel.grid = element_blank()) +
    theme(legend.key.width= unit(8, 'points'), legend.key.height = unit(15, 'points'))
  
  # Extract the legends
  p_legends <- plot_grid(
    NULL,
    cowplot::get_legend(p1 + theme(legend.title = element_text(size = 10), legend.box.margin = margin(t=0, r=0, b=0, l=0))),
    # cowplot::get_legend(p2 + theme(legend.title = element_text(size = 10), legend.box.margin = margin(t=0, r=0, b=0, l=0))),
    cowplot::get_legend(p3 + theme(legend.title = element_text(size = 10), legend.box.margin = margin(t=0, r=0, b=0, l=0))),
    cowplot::get_legend(p4 + theme(legend.title = element_text(size = 10), legend.box.margin = margin(t=0, r=0, b=0, l=0))),
    NULL,
    ncol = 1,
    align = 'v',
    rel_heights = c(1,1.5,2,2,2,1)
  )
  
  p <- plot_grid(p1 + theme(legend.position = "none"), 
                 # p2 + theme(legend.position = "none"), 
                 p3 + theme(legend.position = "none"), 
                 p4 + theme(legend.position = "none"), 
                 p_legends,
                 rel_widths = c(w_p1,1,1,w_legends), 
                 nrow = 1,
                 align = 'h', axis = 'tb')
  
  print(p)
  return(tmp %>% 
           select(dataset, feature_orig, pretty_name, fdr_utest, 
                  increased_in, single_omic_fdr, multi_omic_fdr, 
                  feature_rep_single_omic, feature_rep_multi_omic))
}

# Same style as basic_stats plot, this time describing the diablo components in each dataset
plot_module_stats <- function(sens_analysis_components, 
                              module_cross_view_corrs, 
                              modules_aucs, 
                              feature_type_color_map,
                              dataset_order,
                              hide_y_axis_text = FALSE) {
  
  # ---- Strip 1: N features per dataset per module ---->
  tmp1 <- sens_analysis_components %>%
    mutate(feature_type = substr(feature,1,1)) %>%
    mutate(module = as.numeric(gsub('comp', '', component))) %>%
    group_by(dataset, module, feature_type) %>%
    summarise(N = n(), .groups = "drop") 
  
  tmp1$module2 <- factor(paste('Module', tmp1$module),
                            levels = paste('Module', max(tmp1$module):1))
  
  p1 <- ggplot(tmp1 %>%
                 filter(dataset %in% dataset_order) %>%
                 mutate(dataset = factor(dataset, levels = dataset_order)), 
               aes(x = module2, y = N, fill = feature_type)) +
    geom_bar(color="black", alpha = 0.2, width = 0.8, 
             position="stack", stat="identity") + # Plotted twice as patch
    geom_rect(data = tmp1 %>% filter(dataset %in% dataset_order) %>% 
                mutate(dataset = factor(dataset, levels = dataset_order)) %>%
                select(dataset) %>% distinct(),
              aes(alpha = dataset %>% as.numeric() %% 2 == 0),
              xmin = -Inf, xmax = Inf, ymin = -Inf,
              ymax = Inf, fill = 'gray80', inherit.aes = FALSE) +
    scale_alpha_manual(values = c('FALSE' = 0, 'TRUE' = 0.3), guide = "none") +
    geom_bar(color="black", alpha = 0.5, width = 0.8, 
             position="stack", stat="identity") +
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
    
  # ---- Strip 2: AUC ---->
  tmp2 <- modules_aucs %>%
    filter(run == 'true') %>%
    left_join(
      modules_aucs %>%
        filter(run != 'true') %>%
        group_by(dataset, module) %>%
        summarise(N = n(), 
                  mean_shuf_auc = mean(`auc`, na.rm = TRUE),
                  sd_shuf_auc = stats::sd(`auc`, na.rm = TRUE),
                  .groups = 'drop'),
      by = c('dataset','module')
    ) %>%
    mutate(module = as.numeric(gsub('comp', '', module))) 
  
  tmp2$module2 <- factor(paste('Module', tmp2$module),
                            levels = levels(tmp1$module2))
  points_size <- ifelse(hide_y_axis_text, 3, 3.7)
  
  p2 <- ggplot(tmp2 %>%
                 filter(dataset %in% dataset_order) %>%
                 mutate(dataset = factor(dataset, levels = dataset_order)), 
               aes(x = module2)) +
    geom_rect(data = tmp1 %>% filter(dataset %in% dataset_order) %>% 
                mutate(dataset = factor(dataset, levels = dataset_order)) %>%
                select(dataset) %>% distinct(), 
              aes(alpha = dataset %>% as.numeric() %% 2 == 0), 
              xmin = -Inf, xmax = Inf, ymin = -Inf, 
              ymax = Inf, fill = 'gray80', inherit.aes = FALSE) +
    scale_alpha_manual(values = c('FALSE' = 0, 'TRUE' = 0.3), guide = "none") +
    geom_hline(yintercept = 0.5, color = "darkred", linetype = "dashed", linewidth = 1) +
    # Shuffled standard deviations
    geom_linerange(aes(ymax = mean_shuf_auc + sd_shuf_auc, 
                       ymin = mean_shuf_auc - sd_shuf_auc), 
                   alpha = 0.4, linewidth = 2, color = "grey70") +
    # Shuffled AUCs
    geom_point(aes(y = mean_shuf_auc), 
               shape = 16, size = points_size - 0.5, 
               color = 'grey60', alpha = 0.8) +
    # True standard deviations
    # geom_linerange(aes(ymax = sdev_high_, ymin = sdev_low_, color = feature_set_type),
    #                position = position_dodge(width = p3_dodging), 
    #                alpha = 0.35, linewidth = 2) + 
    # True AUCs
    geom_point(aes(y = `auc`), fill = 'skyblue4',
               shape = 23, color = "black", 
               size = points_size, alpha = 0.9) +
    scale_y_continuous(breaks = seq(0.5,1,0.1)) +
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
  
  # ---- Strip 3: Cross-view correlations ---->
  tmp3 <- module_cross_view_corrs %>%
    filter(run == 'true') %>%
    left_join(
      module_cross_view_corrs %>%
        filter(run != 'true') %>%
        group_by(dataset, module) %>%
        summarise(N = n(), 
                  mean_shuf_corr = mean(avg_pears_corr, na.rm = TRUE),
                  sd_shuf_corr = stats::sd(avg_pears_corr, na.rm = TRUE),
                  .groups = 'drop'),
      by = c('dataset','module')
    ) %>%
    mutate(module = as.numeric(gsub('comp', '', module))) 
  
  tmp3$module2 <- factor(paste('Module', tmp3$module),
                            levels = levels(tmp1$module2))
  
  p3 <- ggplot(tmp3 %>%
                 filter(dataset %in% dataset_order) %>%
                 mutate(dataset = factor(dataset, levels = dataset_order)), 
               aes(x = module2)) + 
    geom_rect(data = tmp1 %>% filter(dataset %in% dataset_order) %>% 
                mutate(dataset = factor(dataset, levels = dataset_order)) %>%
                select(dataset) %>% distinct(), 
              aes(alpha = dataset %>% as.numeric() %% 2 == 0), 
              xmin = -Inf, xmax = Inf, ymin = -Inf, 
              ymax = Inf, fill = 'gray80', inherit.aes = FALSE) +
    scale_alpha_manual(values = c('FALSE' = 0, 'TRUE' = 0.3), guide = "none") +
    # Shuffled standard deviations
    geom_linerange(aes(ymax = mean_shuf_corr + sd_shuf_corr, 
                       ymin = mean_shuf_corr - sd_shuf_corr), 
                   alpha = 0.4, linewidth = 2, color = "grey70") +
    # Shuffled correlations
    geom_point(aes(y = mean_shuf_corr), 
               shape = 16, size = points_size - 0.5, color = 'grey60', alpha = 0.8) +
    # True standard deviations
    # geom_linerange(aes(ymax = sdev_high_, ymin = sdev_low_, color = feature_set_type),
    #                position = position_dodge(width = p3_dodging), 
    #                alpha = 0.35, linewidth = 2) + 
    # True correlations
    geom_point(aes(y = avg_pears_corr), fill = 'sienna4',
               shape = 23, color = "black", 
               size = points_size, alpha = 0.9) +
    scale_x_discrete(expand = c(0, 0.5)) +
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
  tmp_rel_widths <- c(9, 3, 3)
  if (hide_y_axis_text) tmp_rel_widths <- c(7.2, 3, 3)
  print(plot_grid(p1, p3, p2, 
                  nrow = 1, 
                  rel_widths = tmp_rel_widths, 
                  align = 'h', axis = 'tb'))
  
  df_summary <- tmp1 %>%
    group_by(dataset, module) %>%
    summarise(n_features = sum(N), .groups = 'drop') %>%
    left_join(tmp2 %>% select(dataset, module, auc, mean_shuf_auc, sd_shuf_auc),
              by = c('dataset', 'module')) %>%
    left_join(tmp3 %>% select(dataset, module, avg_pears_corr, mean_shuf_corr, sd_shuf_corr),
              by = c('dataset', 'module'))
  
  return(df_summary)
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
      rename(mean_auc = mean_comp_rf_auc, sd_auc = sd_comp_rf_auc) %>%
      mutate(label_n_features = paste0(n_features, ' [', n_modules, ']')) %>%
      select(dataset, mean_auc, sd_auc, label_n_features) %>%
      mutate(model = 'rf') %>%
      mutate(pipeline = caption_minttea) %>%
      mutate(shuf = FALSE),
    summary_aucs %>%
      filter(run != 'true') %>%
      group_by(dataset) %>%
      summarize(mean_auc = mean(mean_comp_rf_auc), sd_auc = sqrt(wtd.var(mean_comp_rf_auc^2))) %>%
      mutate(model = 'rf') %>%
      mutate(pipeline = caption_minttea_shuf) %>%
      mutate(shuf = TRUE),
    summary_aucs %>%
      filter(run == 'true') %>%
      rename(mean_auc = mean_comp_glm_auc, sd_auc = sd_comp_glm_auc) %>%
      mutate(label_n_features = paste0(n_features, ' [', n_modules, ']')) %>%
      select(dataset, mean_auc, sd_auc, label_n_features) %>%
      mutate(model = 'logit') %>%
      mutate(pipeline = caption_minttea) %>%
      mutate(shuf = FALSE),
    summary_aucs %>%
      filter(run != 'true') %>%
      group_by(dataset) %>%
      summarize(mean_auc = mean(mean_comp_glm_auc), sd_auc = sqrt(wtd.var(mean_comp_glm_auc^2))) %>%
      mutate(model = 'logit') %>%
      mutate(pipeline = caption_minttea_shuf) %>%
      mutate(shuf = TRUE)
  ) %>%
    filter(dataset %in% datasets_to_focus_on) %>%
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
  
  p_dodging = 0.55
  
  p <- ggplot(tmp %>% filter(model == 'logit'), 
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
    scale_y_continuous(breaks = seq(0.5,1,0.1)) +
    scale_x_discrete(expand = c(0, 0.5)) +
    scale_fill_manual(name = "Pipeline", values = coloring) +
    guides(fill = guide_legend(reverse=TRUE)) +
    coord_flip() +
    theme_classic() +
    xlab(NULL) +
    ylab("Logistic Regression AUC") +
    theme(panel.grid.major.x = 
            element_line(linewidth = 0.5, color = "grey93")) +
    theme(axis.title.x = element_text(size = 11))
  print(p)
}

plot_module_sensitivity_analysis <- function(sens_analysis_modules, base_module, selected_settings, settings_order = NULL) {
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
  
  if(!is.null(settings_order)) df <- df %>% mutate(run_id = factor(run_id, levels = settings_order))
  
  # Organize features order
  feats_order <- sort(base_features) 
  
  # Hack for nice ordering of features
  if(!is.null(settings_order)) {
    for (i in 2:length(settings_order)) {
      tmp_feats <- df %>% 
        filter(run_id == settings_order[i]) %>%
        pull(feature)
      feats_order <- c(feats_order, tmp_feats[! tmp_feats %in% feats_order])
    }
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
  
  p1 <- ggplot(df, aes(y = feature, x = run_id, fill = module_color)) +
    geom_point(size = 2.5, shape = 21, color = 'black') +
    annotate("rect", xmin = 0.6, xmax = 1.6, ymin = 0, ymax = 1+n_distinct(df$feature), alpha = .1, fill = "cadetblue4") +
    scale_x_discrete(position = 'top') +
    scale_fill_manual(values = c('1' = '#9AD2CB', '2' = '#E0D68A', '3' = '#8E443D', '4' = '#511730')) +
    theme_classic() +
    xlab(NULL) +
    ylab(NULL) +
    theme(legend.position = 'none') +
    theme(axis.text.x.top = element_text(angle = 60, hjust = 0)) +
    theme(axis.text.y = element_text(size = 6)) +
    theme(plot.margin = unit(c(10,95,1,1), 'points')) 
  
  p2 <- ggplot(module_color %>%
                 tidyr::complete(run_id, module_color) %>%
                 mutate(AUC = ifelse(is.na(AUC), 0.501, AUC)), 
               aes(y = AUC, x = run_id, fill = module_color, group = module_color)) +
    geom_col(width = 0.7, color = 'black', position = 'dodge') +
    annotate("rect", xmin = 0.5, xmax = 1.5, ymin = 0, ymax = 1, alpha = .1, fill = "cadetblue4") +
    coord_cartesian(ylim = c(0.5, 1), expand = F) +
    scale_fill_manual(values = c('1' = '#9AD2CB', '2' = '#E0D68A', '3' = '#8E443D', '4' = '#511730')) +
    theme_classic() +
    ylab('Module\'s\n1st PC AUC') +
    xlab(NULL) +
    theme(legend.position = 'none') +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    theme(axis.title.y = element_text(size = 10)) +
    theme(plot.margin = unit(c(10,95,1,1), 'points')) 
  
  return(list(p1 = p1, p2 = p2))
}


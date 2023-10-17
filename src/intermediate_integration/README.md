This folder includes all code related to the intermediate integration (__MintTea__) pipeline.

<img src="https://github.com/borenstein-lab/multi_view_integration_analysis/blob/main/docs/wiki_figure.png" width="700">
     
## Instructions - Running MintTea (`MintTea.R`) on your own data

1. Open an R script from which the MintTea function will be executed.

2. Organize your input data in a *single* data.frame object, following these guidelines:
   * Rows represent samples and columns are features;  
   * The dataframe should include two special columns: a column holding sample identifiers and a column holding study groups ("healthy" and "disease" labels);  
   * Features from each view should start with the following prefixes: 'T_' for taxonomy, 'G_' for genes, 'P_' for pathways, and optionally 'M_' for metabolites;  
   * Features in each view should be pre-processed in advance, according to common practices;  
   * It is highly recommended **to remove rare features, and cluster highly correlated features**;  

3. Run the pipeline by cloning the `multi_view_integration_analysis` repository to your local computer, and then sourcing `MintTea.R` into your script. The function `MintTea(data, ...)` can then be called. 
 
4. Optionally, edit the default pipeline parameters. MintTea supports running the pipeline with multiple parameter combinations, to encourage sensitivity analysis and enable the user to check which settings generate the most informative modules. The full list of MintTea paramaters is given below:

     | Parameter            | Description                                                                                                    |
     | -------------------- | -------------------------------------------------------------------------------------------------------------- |
     | `proc_data`          | A single table containing all features of all views. Samples are rows and features are columns. Two special columns expected to be included in the table are a column holding sample identifiers and a column holding study groups ("healthy" and "disease"). Features from each view should start with the a prefix indicating which view are they related to, followed by two underscores, e.g. 'T_' for taxonomy, 'P_' for pathways, and 'M_' for metabolites. Features in each view should be pre-processed according to common practices. It is advised to remove rare features, and cluster highly correlated features. |
     | `study_group_column` | Name of column holding study groups. |
     | `sample_id_column`   | Name of column holding sample identifiers. |
     | `view_prefixes`      | Feature prefixes differentiating the different views, e.g. `c('T','P','M')` (for Taxonomy, Pathways, Metabolites, respectively). All feature names should start with one of the given prefixes followed by two underscores. |
     | `param_diablo_keepX` | Number of features to select from each view, serving as a constraint for the sparse CCA. Note: these are sparsity constraints for the CCA modules, not the final consensus modules. Higher values will produce larger modules. More than one value can be provided if sensitivity analysis is desired. |
     | `param_diablo_design` | A prior on expected relations between different views. Supports values between 0 and 1 (inclusive). 0 indicates no association between views is expected, and modules should maximize association with disease only. 1 indicates expected inter-view association and modules should therefore maximize both disease-association and between-view associations. More than one value can be provided if sensitivity analysis is desired. |
     | `param_n_repeats`     | Number of sCCA repeats on data subsamples. More than one value can be provided if sensitivity analysis is desired. |
     | `param_n_folds`       | Number of folds used for subsetting the data before running sCCA (DIABLO). A value of 10, for example, will divide the data into 10 subsets, and then run CCA on 9/10 of the data, excluding each subset one at a time. Lower values will result in smaller subsets used for training and accordingly to higher variability between sCCA models. In such a case we expect less modules to be identified, but their robustness to be higher. More than one value can be provided if sensitivity analysis is desired. |
     | `param_diablo_ncomp`  | Number of sCCA components to extract each DIABLO run. Note that DIABLO authors recommend using only the first few components. Typically, components >3 are less stable, and will often not contribute to final consensus modules. More than one value can be provided if sensitivity analysis is desired. |
     | `param_edge_thresholds` | Number between 0 and 1 (exclusive), determining the threshold for consensus components. Higher values mean more conservative results. Values between 0.5 and 0.8 are recommended. More than one value can be provided if sensitivity analysis is desired.  |
     | `n_evaluation_repeats` | Number of cross-validation repeats for overall AUROC estimation. |
     | `n_evaluation_folds`  | Number of cross-validation folds for overall AUROC estimation. |
     | `log_level`           | See `library(logger); ?log_levels` |
     | `seed`                | For result replicability. |


6. Pipeline results are returned as a list of multi-view modules, given for each MintTea pipeline setting requested. For each module, the following properties are returned: 

     | Module property           | Details                                                                                                             |
     | ------------------------- | ------------------------------------------------------------------------------------------------------------------- |
     | `module_size`             | The number of features in this module.                                                                              |
     | `features`                | 1st prinicipal component (PC) of each module, for each pipeline setting.                                            |
     | `module_edges`            | Edge weights for every pair of features in this module that co-occured in sGCCA components at least once. Edge weights are calculated as the number of times each pair co-occured in the same sGCCA component, divided by `param_n_repeats` * `param_n_folds`. These weights are given in case the user wants to draw the module as a network.  |
     | `auroc`                   | AUROC of each module by itself, describing the module's association with the disease. Computed using its first PC and evaluated over repeated cross-validation. Note: It is warmly advised to further evaluate module-disease associations using an independent test set. |
     | `shuffled_auroc`          | As above, but using 99 randomly sampled modules of the same size and same proprtions of views.                      |
     | `inter_view_corr`         | Average correlation between features from different views.                                                          |
     | `shuffled_inter_view_corr` | As above, but using 99 randomly sampled modules of the same size and same proprtions of views.                     |
  
13. To evaluate the obtained results, we recommend starting by examining the following:

   * For each pipeline setting - how many modules were found, and what are the module sizes (i.e., number of features included)? 
   * What was the AUC achieved by each module alone? (see `auroc`)
   * How does this AUC compare to the random-modules AUC's?


Tips:
   
   * Optimal module sizes depend on the downstream analysis. For manual interpretation for example, smaller modules may be favorable. If your modules came out too large, consider decreasing `param_diablo_keepX`, or decreasing `param_n_folds`, or increasing `param_edge_thresholds`. Symmetrically, if your modules are too small consider the opposite.
   * If the overall AUC is low, and/or all individual module AUC's are low, you may want to consider decreasing `param_diablo_design`, effectively assigning a higher importance to associations with disease as opposed to associations in-between views.
  
For questions about the pipeline, please contact elbo@tauex.tau.ac.il.

***

**Backlog:**

     * Support parallel running to shorten runtimes.
     * Generalize to support any categorical label (currently hard-coded to 'healthy' and 'disease').
     * Generalize to support continuous labels.
     

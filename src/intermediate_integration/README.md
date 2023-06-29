This folder includes all code related to the intermediate integration (__MintTea__) pipeline.

<img src="[markdownmonstericon.png](https://github.com/borenstein-lab/multi_view_integration_analysis/blob/main/docs/wiki_figure.png)"
     alt="MintTea input and output"
     style="float: left; margin-right: 10px;" />
     
## Instructions - Running MintTea (`MintTea.R`) on your own data

1. Open an R script where MintTea function will be executed.

2. Organize your input data in a *single* data.frame object, following these guidelines:
   * Rows represent samples and columns are features;  
   * The dataframe should include two special columns: a column holding sample identifiers and a column holding study groups ("healthy" and "disease" labels);  
   * Features from each view should start with the following prefixes: 'T_' for taxonomy, 'G_' for genes, 'P_' for pathways, and optionally 'M_' for metabolites;  
   * Features in each view should be pre-processed in advance, according to common practices;  
   * It is advised to remove rare features, and cluster highly correlated features;  

3. Run the pipeline by copying the `MintTea.R` script locally, sourcing it into your script, and then calling the `MintTea(data, ...)` function. 
 
4. Optionally, edit the default pipeline parameters. MintTea supports running the pipeline with multiple parameter combinations, to encourage sensitivity analysis and enable the user to check which setting generated most informative modules. The full list of MintTea paramaters is given below:

     | Parameter            | Description                                                                                                    |
     | -------------------- | -------------------------------------------------------------------------------------------------------------- |
     | `proc_data`          | A single table containing all features of all views. Samples are rows and features are columns. Two special columns expected to be included in the table are a column holding sample identifiers and a column holding study groups ("healthy" and "disease"). Features from each view should start with the following prefixes: 'T_' for taxonomy, 'G_' for genes, 'P_' for pathways, and optionally 'M_' for metabolites (future versions will support arbitrary prefixes). Features in each view should be pre-processed according to common practices. It is advised to remove rare features, and cluster highly correlated features. |
     | `study_group_column` | Name of column holding study groups. |
     | `sample_id_column`   | Name of column holding sample identifiers. |
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


6. Pipeline results are all returned in a single R list, and contain both the detailed modules and various evaluations and statistics about them. Main outputs include:

     | Output                    | Description                                                                                                         |
     | ------------------------- | ------------------------------------------------------------------------------------------------------------------- |
     | `sens_analysis_modules`   | The table lists all identified modules (i.e., the full list of features in each module), for each pipeline setting. |
     | `latent_vars`             | 1st prinicipal component (PC) of each module, for each pipeline setting.                                            |
     | `module_variance_expl`    | Variance explained by first PC of the true modules, as well as shuffled modules. True modules are expected to have significantly higher levels of variance explained in their first PC (compared to shuffled modules) as features are highly associated with one another. |
     | `module_inter_omic_cors`  | Average correlation between features from different views, per module and per pipeline setting. Again, compared to shuffled modules. |
     | `summary_module_aucs`     | AUROC of each module by itself, describing the module's association with the disease. Computed using its first PC and evaluated over repeated cross-validation. Compared to shuffled modules. |
     | `summary_overall_aucs`    | Combined AUROC of all modules of a dataset, using their first PC and a simple random forest or logistic regression model, and evaluated over repeated cross-validation. Describes the overall predictive power of all modules combined. Compared to shuffled modules. |
  
7. To evaluate the obtained results, we recommend examining the following:

   * XXX
   * XXX
   * XXX
  
For questions about the pipeline, please contact elbo@tauex.tau.ac.il.

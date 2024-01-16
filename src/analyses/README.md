This folder includes all code related to the analysis and visualization of the intermediate integration pipeline outputs, as presented in the manuscript.

Specifically, it includes the following scripts/notebooks: 

_Analysis of Random forest early integration results_

* `rf_results_analysis.Rmd`: Notebook to summarize ML results on each dataset.

_Comparisons between MintTea and other sCCA methods_

* `false_discovery_analysis.R`: Compares MintTea to MultiCCA and DIABLO in terms of their false discovery rates, i.e. their tendency to detect valid modules even in randomly shuffled data.  
* `minttea_comparisons_module_associations.R`: Compares MintTea to MultiCCA and DIABLO in terms of the predictivity of each individual module/component and in terms of the correlation between features from different views.  
* `minttea_comparisons_stability_analysis.R`: Compares MintTea to MultiCCA and DIABLO in terms of stability of the obtained components/modules.  
* `mintea_comparisons.Rmd`: A notebook that visualizes and summarizes all the various comparisons performed above, using the outputs of the above scripts.   

_Main analysis of MintTea results on multiple datasets_

* `intermediate_integration_analysis.Rmd`: Main analysis notebook. Includes most of the manuscript figures and tables.

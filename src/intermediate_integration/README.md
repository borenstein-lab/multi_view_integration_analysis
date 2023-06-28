This folder includes all code related to the intermediate integration (__MintTea__) pipeline.

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

     | Parameter | Description |
     | --------- | ----------- |
     | XXX       | YYY         |

5. Pipeline results are all returned in a single R list, and contain the following outputs:

     | Output    | Description |
     | --------- | ----------- |
     | XXX       | YYY         |

6. To evaluate the obtained results, we recommend examining the following:

   * XXX
   * XXX
   * XXX
  
For questions about the pipeline, please contact elbo@tauex.tau.ac.il.

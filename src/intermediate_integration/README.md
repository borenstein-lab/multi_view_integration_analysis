This folder includes all code related to the intermediate integration pipeline.

## Instructions - Running the pipeline on a single dataset

1. Prepare input data as follows:

1.1. metadata file with 'sample_id' as the first column and 'DiseaseState' as an additional column in which one of two labels (e.g. 'Control' & 'Disease') appear per sample.
1.2. taxonomy file with 'sample_id' as the first column and all other columns representing different taxa and their abundances (could be relative abundances or counts).
1.3. pathways file with 'pathway' as the first column (pathway names/codes), and other columns corresponding to the sample_ids.
1.4. gene family files with 'function' as the first column (gene family name/code), and other columns corresponding to the sample_ids.
1.5. (optional) metabolite file with 'sample_id' as the first column and all other columns representing different metabolites and their abundances.

All 4 (or 5) files should be placed in one folder. The folder's name is considered as the dataset's name. File names are not restricted, they just need to be updated in the configuration file (see next step).

2. Go over the `config.yml` file an update all fields as needed. See documentation within the config file for a description of each field. Make sure that the paths to the input files prepared in the previous step are correct.

3. Edit the names of output files/folders as needed, under the "Define output files" title in the `intermediate_integration_main.R` script.

4. Optionally, edit the pipeline parameters you want to test in the lines after the title "Run sensitivity analysis on: ...". The pipeline by default uses a few different settings to encourage sensitivity analysis and enable the user to check which setting generated most informative modules.

5. Run the pipeline on a single dataset either by using command line (see below) or by running the script from within R studio.
`Rscript intermediate_integration_main.R -d <dataset name> -w <working directory>`

The pipeline uses all of the following scripts:
`diablo_utils.R` - Misc. utilities
`ml_pipeline/preprocessing.R` - Data pre-processing functions
`ml_pipeline/clustering.R` - Clustering of highly redundant features

6. Outputs are all given in a single RData file. See `intermediate_integration_analysis.Rmd` notebook for optional follow-up analysis.


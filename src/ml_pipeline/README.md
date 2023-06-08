This folder includes all code related to the machine learning (ML) pipeline.

## Instructions - Running the pipeline on a single dataset

1. Prepare input data as follows:

1.1. metadata file with 'sample_id' as the first column and 'DiseaseState' as an additional column in which one of two labels (e.g. 'Control' & 'Disease') appear per sample.
1.2. taxonomy file with 'sample_id' as the first column and all other columns representing different taxa and their abundances (could be relative abundances or counts).
1.3. pathways file with 'pathway' as the first column (pathway names/codes), and other columns corresponding to the sample_ids.
1.4. gene family files with 'function' as the first column (gene family name/code), and other columns corresponding to the sample_ids.
1.5. (optional) metabolite file with 'sample_id' as the first column and all other columns representing different metabolites and their abundances.

All 4 (or 5) files should be placed in one folder. The folder's name is considered as the dataset's name. File names are not restricted, they just need to be updated in the configuration file (see next step).

2. Go over the `config.yml` file an update all fields as needed. See documentation within the config file for a description of each field. Make sure that the paths to the input files prepared in the previous step are correct.

3. Run the pipeline on a single dataset either by using command line (see below) or by running the script from within R studio.
`Rscript ml_pipeline.R -d <dataset name> -w <working directory>`

Note that the pipeline can also be run in a "test" mode using the "test" part of the configuration file (`config.yml`)

The pipeline uses all of the following scripts:
`utils.R` - Misc. utilities
`preprocessing.R` - Data pre-processing functions
`feature_selection.R` - Feature selection methods
`clustering.R` - Clustering of highly redundant features
`postprocess_results.R` - Functions for post processing of ML results

4. Lastly, see `postprocess_results.R` for various optional functions for summarizing ML results.


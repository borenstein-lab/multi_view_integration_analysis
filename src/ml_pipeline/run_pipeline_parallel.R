library(config)

source('utils.R')

run_pipeline_in_parallel <- function(dataset_names = NULL) {
  
  # Create all output directories
  utils_create_output_dirs()
  
  # Create template for command 
  script_cmd_template <- "Rscript ml_pipeline.R -d %s"
  
  # If a list of datasets is provided, take it. Otherwise, use all datasets in 'ml_input'
  if (is.null(dataset_names)) {
    dataset_names <- utils_get_available_ds()
  } 
  
  # Run the ml_pipeline script per dataset
  results <- lapply(dataset_names, function(ds_name) {
    log_path <- sprintf(config::get('paths_templates')$log, ds_name)
    script_cmd <- sprintf(script_cmd_template, ds_name)
    script_cmd <- paste(script_cmd, "2>", log_path)
    print(script_cmd)
    
    return(system(script_cmd,
                  wait = FALSE,
                  show.output.on.console = FALSE))
  })
}

###################################
# Main 
###################################

# Run pipeline on all datasets, by default - all datasets found in the ml_input package
run_pipeline_in_parallel()

# Tests
# run_pipeline_in_parallel(dataset_names = c('t1d_alkanani'))


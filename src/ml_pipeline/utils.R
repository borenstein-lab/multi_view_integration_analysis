require(config)
require(logger)
require(readr)

options(stringsAsFactors = FALSE)  # For compatibility with R < 4


utils_write_excel <- function(df, row.names = FALSE, col.names = TRUE, ...) {
  write.table(df,
              "clipboard",
              sep = "\t",
              row.names = row.names,
              col.names = col.names,
              ...)
}


utils_create_output_dirs <- function() {
  dir.create(config::get('paths')$ml_output_dir, showWarnings = FALSE)
  dir.create(config::get('paths')$logs_dir, showWarnings = FALSE)
  dir.create(config::get('paths')$results_tables_dir, showWarnings = FALSE)
}


utils_save_tsv_table <- function(data, path, quote = FALSE, sep = '\t', append = FALSE, col.names = TRUE) {
  write.table(data, file = path, quote = quote, sep = sep, append = append, row.names = FALSE, col.names = col.names)
}


utils_get_available_ds <- function() {
  return(list.dirs(path = config::get('paths')$ml_input_dir,
                   recursive = FALSE,
                   full.names = FALSE))
}

utils_get_kegg_compound_names_mapping <- function() {
  # Read mapping from KEGG compound ID's to names
  kegg_metab_names <- read_delim(config::get("paths")$kegg_metabolites_names, 
                                 delim = "\t", 
                                 escape_double = FALSE, 
                                 col_names = c("id","name"), 
                                 show_col_types = FALSE,
                                 trim_ws = TRUE)
  
  # If there are several names, take the first
  kegg_metab_names <- kegg_metab_names %>%
    rename(all_names = name) %>%
    tidyr::separate(all_names, 
             into = c("name"), 
             sep = ";", 
             remove = FALSE, 
             extra = "drop")
  
  ## Convert into a map (named vector instead of a table)
  kegg_metab_names_map <- kegg_metab_names$name
  names(kegg_metab_names_map) <- kegg_metab_names$id
  
  return(kegg_metab_names_map)
}

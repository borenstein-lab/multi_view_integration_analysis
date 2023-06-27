# Pre-processing of the pathways DB from MetaCyc.
# After running this file, an object called "metacyc_db" is available for querying the pathways database.
# Additionally, a list of bacteria-associated pathways can be used to filter non-relevant pathways (common artefacts)

library(config)
library(logger)
library(dplyr)
library(collections)
library(readr)

# Read configuration file
config <- config::get(file = "src/ml_pipeline/config.yml")

# Reads the database from MetaCyc database path. The database is a text file
# generated from MetaCyc smart tables (available in the website after registration).
# The database file should contain at least two columns: "Pathways", "Super-Pathways".
# The function generates three data structures:
# 1. leaves - a set of leaves in the database, i.e. the lower pathways in the pathways hierarchy.
# 2. sub_to_super - a dictionary with a mapping between a subpathway to its superpathways set.
# 3. super_to_sub - a dictionary with a mapping between a superpathway to its subpathways set.
read_database <- function(metacyc_db_path) {
  data <- read.csv(file = metacyc_db_path,
                   sep = '\t', 
                   header = TRUE, 
                   quote = "", 
                   check.names = FALSE)
  
  # Build subpathway -> superpathways mapping
  all_superpathways <- sets::set()
  
  sub_to_super <- dict()
  have_superpathway <- data %>%
    filter(`Super-Pathways` != '')
  
  for (i in 1:nrow(have_superpathway)) {
    subpathway <- have_superpathway[i, 'Pathways']
    superpathways_raw <- have_superpathway[i, 'Super-Pathways']
    superpathways <- strsplit(superpathways_raw, ' // ')[[1]] %>% sets::as.set()
    
    sub_to_super$set(subpathway, superpathways)
    all_superpathways <- sets::set_union(all_superpathways, superpathways)
  }
  
  # Build leaves set
  leaves <- data %>%
    select(`Pathways`) %>%
    unlist() %>%
    unname() %>%
    sets::as.set() %>%
    sets::set_symdiff(all_superpathways)
  
  # Build superpathway -> subpathways mapping
  super_to_sub <- dict()
  for (subpathway_name in names(sub_to_super$as_list())) {
    for (superpathway_name in sub_to_super$get(subpathway_name)) {
      if (!super_to_sub$has(superpathway_name)) {
        super_to_sub$set(superpathway_name, sets::set())
      }
      super_to_sub$set(superpathway_name,
                          sets::set_union(super_to_sub$get(superpathway_name),
                                          sets::set(subpathway_name)))
    }
  }
  
  return(list(leaves=leaves,
              sub_to_super=sub_to_super,
              super_to_sub=super_to_sub))
}


# Receives a single pathway and a set of the pathways in the data frame and decides
# if the pathway should filtered out or stay. The reason for filtering pathways
# is that there is a pathway hierarchy in MetaCyc, leading to both pathway and
# it's sub-pathways to exist in the dataset's pathways. This condition is problematic
# since converting to relative abundance will result in false abundances due to
# the duplicate pathways. The general logic is to keep only the lower level of
# pathways (no super-pathways), the exact logic is documented inside the function.

should_keep_pathway <- function(ds_pathway, all_df_pathways) {
  # If it's a leaf in MetaCyc path hierachy (doesn't have sub pathways), keep it
  if (sets::set_contains_element(metacyc_db$leaves, ds_pathway)) {
    return(TRUE)
  }
  
  # If MetaCyc doesn't contain this pathway, keep it, probably the database isn't updated. 
  if (!metacyc_db$super_to_sub$has(ds_pathway)) {
    return(TRUE)
  }
  
  # We already know it's not a leaf in MetaCyc heirachy, but it can be a leaf
  # in this particular dataset, because the dataset doesn't necessarily contain
  # all the pathways. So now check if it's a superpathway without children in
  # this dataset, practically a leaf in relation to this dataset.
  children <- metacyc_db$super_to_sub$get(ds_pathway)
  if (length(sets::set_intersection(children, all_df_pathways)) == 0) {
    return(TRUE)
  }
  return(FALSE)  
}

# Called from pre-processing.R, receives a dataframe with a column named "pathway"
# and filters the pathways according to the logic in "should_keep_pathway" function.
filter_super_pathways <- function(pathways_df) {
  ds_pathways <- pathways_df %>%
    select(pathway) %>%
    unlist() %>%
    unname() %>%
    sets::as.set()
  
  ds_pathways_list <- as.list(ds_pathways)
  names(ds_pathways_list) <- ds_pathways_list
  
  results <- sapply(ds_pathways_list, 
                    should_keep_pathway,
                    all_df_pathways = ds_pathways)
  
  filtered_df <- pathways_df %>%
    filter(pathway %in% names(results[results]))
  
  log_debug(sprintf('Pathways DB: filtered %d/%d redundant superpathways',
                  (length(ds_pathways) - nrow(filtered_df)),
                  length(ds_pathways)))
  return(filtered_df)
}

# Return a list of bacterial MetaCyc pathways only
get_bacterial_pathways <- function() {
  require(readr)
  
  all_pwys <- read_delim(config$paths$metacyc_pathways2taxa, 
                         delim = "\t", escape_double = FALSE, 
                         col_select = c('Object ID', 'Taxonomic-Range'),
                         show_col_types = FALSE, trim_ws = TRUE)
  all_pwys <- all_pwys %>%
    rename(pwy = 1, tax = 2) %>%
    tidyr::separate(col = tax, into = paste0('tax', 1:30), sep = ' // ', fill = 'right') %>%
    tidyr::pivot_longer(cols = -pwy, names_to = 'dummy', values_to = 'tax', values_drop_na = TRUE) %>%
    select(-dummy)
  
  all_bacteria <- read_delim(config$paths$metacyc_bacterial_taxa, 
                             delim = "\t", escape_double = FALSE, 
                             col_select = c('Object ID', 'Common-Name'), 
                             show_col_types = FALSE, trim_ws = TRUE)
  all_bacteria <- all_bacteria %>%
    rename(tax = 1, tax_name = 2)
  
  # Filter metacyc pathways to those associated with bacteria only
  all_pwys <- all_pwys %>%
    filter(tax %in% all_bacteria$tax)
  
  return(unique(all_pwys$pwy))
}

filter_bacterial_pathways <- function(pathways_df) {
  bacterial_pwys <- get_bacterial_pathways()
  filtered_df <- pathways_df %>%
    filter(pathway %in% bacterial_pwys)
  log_debug(sprintf('Pathways DB: filtered %d/%d non-bacterial pathways',
                    (nrow(pathways_df) - nrow(filtered_df)), 
                    nrow(pathways_df)))
  return(filtered_df)
}

# Load the database only once
if (!exists('metacyc_db')) {
  metacyc_db <- read_database(config$paths$metacyc_pathways_hierarchy)
}



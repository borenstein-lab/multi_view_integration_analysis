n.nodes <- 4

# Server settings
#.libPaths("/specific/elhanan/PROJECTS/16S_BETTER_PREDS_IS/R4_env")
# n.nodes <- 50

library(config)
library(dplyr)
library(readr)
library(S4Vectors)
library(fs)
library(Matrix.utils)
library(foreach)
library(doSNOW)
library(curatedMetagenomicData)

options(scipen = 999)

cl <- makeCluster(n.nodes)
registerDoSNOW(cl)

# We search for relevant studies that meet the criteria:
#  (1) At least 40 cases and at least 40 controls
#  (2) Metadata, taxonomy, ko's and pathways data - available

# relevant_study_groups <- sampleMetadata %>%
#   group_by(study_name, study_condition) %>%
#   summarise(N = n(), N_unique_individuals = n_distinct(subject_id)) %>%
#   filter(N_unique_individuals >= 30) %>%
#   group_by(study_name) %>%
#   filter(n() >= 2) %>%
#   ungroup() %>%
#   group_by(study_name) %>%
#   summarise(N_study_groups = n(),
#             N_total_individuals = sum(N_unique_individuals),
#             study_conditions = paste0(study_condition, collapse = ", "))
# 
# # TODO: in some longitudinal studies samples "move" between study conditions. Need to propely deal with these cases, for now exclude these studies.
# longitudinal_studies_with_subject_conversions <- sampleMetadata %>%
#   group_by(study_name, subject_id) %>%
#   summarise(N = n(), N_conditions = n_distinct(study_condition)) %>%
#   filter(N_conditions > 1)


# Given the path to a ko <--> UniRef mapping table (from humann), load the table into R.
# This function leaves the table in a wide format, i.e. the first column represents a ko, 
#  and the following columns represent all Uniref codes associated with it
get_uniref_ko_mapping <- function(mapping_path) {
  uniref_mapping <- read_delim(mapping_path, 
                               delim = "\t", escape_double = FALSE, 
                               show_col_types = FALSE,
                               col_names = c("KO", "UniRef90"), 
                               trim_ws = TRUE)
  
  uniref_mapping$UniRef90_short <- gsub("UniRef90_","",uniref_mapping$UniRef90)
  
  return(uniref_mapping)
}



# Given the original study name and the two study groups to include, the function returns 
#  the relevant subset of sampleMetadata (which is the entire metadata of CuratedMetagenomicData)
# Additionally, when multiple samples per subject exist, only the first is taken
get_samples_to_include <- function(study_name_orig, group1, group2) {
  samples_to_include <- sampleMetadata %>%
    filter(study_name == study_name_orig) %>%
    filter(study_condition %in% c(group1, group2)) %>%
    select(where(~ !all(is.na(.x))))
  
  # Keep only one sample per subject (if possible, take the first one, assuming it represents a baseline sample)
  if ("days_from_first_collection" %in% names(samples_to_include)) {
    samples_to_include <- samples_to_include %>%
      arrange(days_from_first_collection)
  }
  samples_to_include <- samples_to_include %>%
    group_by(subject_id) %>%
    filter(row_number()==1) %>%
    ungroup()
  
  return(samples_to_include)
}



# Prepare taxonomy profiles in expected format
prepare_dataset_taxonomy <- function(samples_to_include, output_path) {
  taxonomy <- returnSamples(samples_to_include, "relative_abundance", rownames = "short")
  taxonomy <- data.frame(t(assay(taxonomy)), check.names = FALSE)
  taxonomy <- tibble::rownames_to_column(taxonomy, "sample_id")
  write.table(taxonomy, file = output_path, row.names = FALSE, sep = "\t", quote = FALSE)
}



# Prepare pathway profiles in expected format
# Note: pathways/genes are given in RPK units
prepare_dataset_pathways <- function(samples_to_include, output_path) {
  pwys <- returnSamples(samples_to_include, "pathway_abundance")
  pwys <- data.frame(assay(pwys), check.names = FALSE)
  
  # Collapse to get un-stratified profiles
  pwys$pathway <- gsub("\\|.*$", "", rownames(pwys))
  pwys <- aggregate(.~pathway, pwys, sum)
  
  # Split pathway name into code+description
  pwys <- pwys %>%
    tidyr::separate(col = pathway, into = c("pathway", "description"), sep = ": ", fill = "right") %>% 
    relocate(pathway, description)
  
  ## Not needed: Transpose in order to get samples in rows and pathways in columns
  #pwys <- pwys %>% tibble::column_to_rownames("unstrat_pathway")
  #pwys <- data.frame(t(pwys), check.names = FALSE)
  #pwys <- tibble::rownames_to_column(pwys, "sample_id")
  
  write.table(pwys, file = output_path, row.names = FALSE, sep = "\t", quote = FALSE)
}



# Prepare KO profiles in expected format (for cmp)
prepare_dataset_genes <- function(dataset_desc, samples_to_include, output_path, uniref_mapping) {
  # 1. Get gene abundances, given in UniRef and including both stratified and unstratified outputs
  # 1.1. If an rda file is available - take it from there
  if (!is.na(dataset_desc$gene_abundances_rda)) {
    loaded <- load(dataset_desc$gene_abundances_rda)
    gene_families <- get(loaded[1])
    gene_families <- gene_families[,samples_to_include$sample_id]
    # Keep only relevant samples
  } else { # 1.2. Else, download from ExperimentHub directly
    gene_families <- returnSamples(samples_to_include, "gene_families")
    gene_families <- assay(gene_families)
  }
  message("  Got gene families")
  
  # 2. Keep only unstratified profiles
  unstratified_rows <- !grepl("\\|", rownames(gene_families))
  gene_families_unstrat <- gene_families[unstratified_rows,]
  rm(gene_families) # Free space
  message("  Kept only unstratified profiles")
  
  # 3. Replace rownames into short UniRef names (delete the "UniRef90_" prefix), for efficiency
  uniref_short_names <- gsub("UniRef90_","",rownames(gene_families_unstrat))
  rownames(gene_families_unstrat) <- uniref_short_names
    
  # 4. Create a compact map from ko's to associated UniRefs 
  #  (i.e. a named list where the name is the ko code and the value is a vector of uniref's)
  compact_uniref_mapping <- uniref_mapping %>% filter(UniRef90_short %in% uniref_short_names) %>% arrange(UniRef90_short)
  ko_to_unirefs <- as.list(compact_uniref_mapping$UniRef90_short)
  names(ko_to_unirefs) <- compact_uniref_mapping$KO
  ko_to_unirefs = tapply(unlist(ko_to_unirefs, use.names = FALSE), 
                         rep(names(ko_to_unirefs), lengths(ko_to_unirefs)), FUN = c)
  message("  Created compact ko-uniref mapping")
  
  # 5. Regroup the table
  # For each ko, sum all relevant rows in current matrix
  # Note: UniRef's associated with more than one ko will be summed more than once, 
  #  mimicing humann's default regrouping algorithm
  pb <- txtProgressBar(style = 3, max = length(ko_to_unirefs)) # Custom progress bar
  opts_snow <- list(progress = function(n) setTxtProgressBar(pb, n))
  gene_families_ko <- foreach(i = 1:length(ko_to_unirefs), 
                              .combine = 'rbind', 
                              .packages=c("Matrix"),
                              .options.snow = opts_snow) %dopar% {
    ko <- names(ko_to_unirefs)[i]
    unirefs_to_group <- ko_to_unirefs[[ko]]
    if (length(unirefs_to_group) == 1) {
      gene_families_unstrat[unirefs_to_group,]
    } else {
      Matrix::colSums(gene_families_unstrat[unirefs_to_group,])
    }
  }
  rownames(gene_families_ko) <- names(ko_to_unirefs)
  message("  Performed the mapping")
  
  # Dump to file
  gene_families_ko_df <- data.frame(gene_families_ko, check.names = FALSE)
  gene_families_ko_df <- tibble::rownames_to_column(gene_families_ko_df, "function")
  write.table(gene_families_ko_df, file = output_path, row.names = FALSE, sep = "\t", quote = FALSE)
}



# Prepare sample metadata in expected format
preparet_dataset_metadata <- function(samples_to_include, output_path) {
  metadata <- samples_to_include %>%
    dplyr::rename(DiseaseState = study_condition) %>%
    mutate(DiseaseState = ifelse(DiseaseState == "control", "H", DiseaseState)) %>%
    relocate(sample_id)
  
  write.table(metadata, file = output_path, row.names = FALSE, sep = "\t", quote = FALSE)
}



# Prepare all files needed for further processing of the dataset
prepare_dataset <- function(dataset_desc, uniref_mapping, override = FALSE) {
  message(sprintf('Preparing dataset %s', dataset_desc$new_study_name))
  
  # Extract original study name and wanted study groups
  study_name_orig <- dataset_desc$study_name
  group1 <- dataset_desc$study_condition_1
  group2 <- dataset_desc$study_condition_2
  dataset_name <- dataset_desc$new_study_name
  
  # Define paths for output files
  metadata_path <- sprintf(config::get('paths_templates')$metadata, dataset_name)
  tax_path <- sprintf(config::get('paths_templates')$taxonomy, dataset_name)
  pwy_path <- sprintf(config::get('paths_templates')$pathways, dataset_name)
  gene_path <- sprintf(config::get('paths_templates')$metagenome, dataset_name)
  
  # Create output directory if do not exist yet (does nothing if they do exist)
  dir.create(path_dir(metadata_path), showWarnings = FALSE)
  
  # Get a subset of metadata including relevant samples
  samples_to_include <- get_samples_to_include(study_name_orig, group1, group2)
  message(sprintf('Dataset includes %i controls and %i cases', 
                  nrow(samples_to_include %>% filter(study_condition == group1)),
                  nrow(samples_to_include %>% filter(study_condition == group2))))
  
  # Prepare metadata in expected format
  if (override | (! file.exists(metadata_path))) preparet_dataset_metadata(samples_to_include, output_path = metadata_path)
  message("Done - metadata")
  
  # Prepare taxonomy profiles in expected format
  if (override | (! file.exists(tax_path))) prepare_dataset_taxonomy(samples_to_include, output_path = tax_path)
  message("Done - taxonomy")
  
  # Prepare pathway abundances in expected format
  if (override | (! file.exists(pwy_path))) prepare_dataset_pathways(samples_to_include, output_path = pwy_path)
  message("Done - pathways")
  
  # Prepare UniRef genes dataset and convert it to KO profiles for CMP calculation
  if (override | (! file.exists(genes_for_cmp_path))) 
    prepare_dataset_genes(dataset_desc, samples_to_include, output_path = gene_path, uniref_mapping = uniref_mapping)
  message("Done - gene families")
  
  message("Done\n")
}


############
#   Main   #
############

# Load the configuration specifying which studies to use and which study groups to take from each study
studies_to_process <- read_csv(config::get('paths')$curatedMetagenomicData_study_list, comment = "#", show_col_types = FALSE)
studies_to_process_l <- split(studies_to_process, 1:nrow(studies_to_process))
names(studies_to_process_l) <- studies_to_process$new_study_name

# Get a mapping table from uniref90 to ko's
uniref_mapping <- get_uniref_ko_mapping(config::get('paths')$humann_uniref_ko_map)

# For testing:
# dataset_name <- "cirrhosis_qin_2014"
# dataset_desc <- studies_to_process_l[[dataset_name]]
# Or:
# studies_to_process_l <- studies_to_process_l[11:16]

# Iterate over all desired datasets, prepare data files for each one
dummy_result <- lapply(studies_to_process_l, prepare_dataset, uniref_mapping = uniref_mapping)

# Clean up
stopCluster(cl)
message("Cluster stopped, script completed")
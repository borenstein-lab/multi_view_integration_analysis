require(config)
require(logger)
require(readr)

utils_create_output_dirs <- function() {
  dir.create(config$paths$ml_output_dir, showWarnings = FALSE)
  dir.create(config$paths$logs_dir, showWarnings = FALSE)
  dir.create(config$paths$results_tables_dir, showWarnings = FALSE)
}


utils_save_tsv_table <- function(data, path, quote = FALSE, sep = '\t', append = FALSE, col.names = TRUE) {
  write.table(data, file = path, quote = quote, sep = sep, append = append, row.names = FALSE, col.names = col.names)
}


utils_get_available_ds <- function() {
  return(list.dirs(path = config$paths$ml_input_dir,
                   recursive = FALSE,
                   full.names = FALSE))
}

# Combine features' importance p values using Fisher's method 
post_combine_p_vals_fisher <- function(pvals) {
  if (length(pvals) == 1) return(pvals)
  if (sum(!is.na(pvals)) == 0) return(NA)
  return(sumlog(p.adjust(pvals, method="fdr"))$p)
}
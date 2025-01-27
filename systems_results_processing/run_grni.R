#!/usr/bin/env Rscript
# run_grni.R

library(gurobi)
library(GRNOpt)
library(stringr)
library(dplyr)
library(readr)

remove_duplicates <- function(df, id_cols) {
  n_before <- nrow(df)
  df_unique <- df %>% distinct(across(all_of(id_cols)), .keep_all = TRUE)
  n_after <- nrow(df_unique)
  if (n_before > n_after) {
    warning(paste("Removed", n_before - n_after, "duplicate rows"))
  }
  return(df_unique)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript run_grni.R <bool_folder> <pkn_path> <output_folder>")
}

bool_folder <- args[1]
pkn_path <- args[2]
output_folder <- args[3]

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

choicerule <- "ID"

pkn <- tryCatch({
  read.csv(pkn_path, sep = ',', header = TRUE, stringsAsFactors = FALSE, quote = "")
}, error = function(e) {
  stop(paste("Error reading PKN file:", e$message))
})

colnames(pkn) <- c("TF", "effect", "target")
pkn <- remove_duplicates(pkn, c("TF", "effect", "target"))

CellTypes <- list.files(bool_folder, pattern = "\\.csv$", full.names = TRUE)

for (bool_file in CellTypes) {
  tryCatch({
    file_name <- tools::file_path_sans_ext(basename(bool_file))
    output_file <- file.path(output_folder, paste0('GRNI_', file_name, '.csv'))
    
    # Read boolean data and remove duplicates
    bool_data <- read.csv(bool_file, sep = ',', header = TRUE, stringsAsFactors = FALSE, quote = "")
    print(head(bool_data))
    colnames(bool_data) <- c("TF", "TF_value")
    bool_data <- remove_duplicates(bool_data, "TF")
    
    duplicates <- bool_data %>% group_by(TF) %>% filter(n() > 1)
    if (nrow(duplicates) > 0) {
      warning(paste("Duplicates found in bool_data after removal:", 
                    paste(duplicates$TF, collapse = ", ")))
    }
    
    bool_tfonly <- bool_data[bool_data$TF %in% unique(pkn$TF) | bool_data$TF %in% unique(pkn$target), ]
    
    result <- prune_gurobi(pkn, bool_tfonly, rule = choicerule)
    
    print(paste("Processing:", file_name))
    print(result)
    
    #check for duplicates in result before writing
    result_unique <- remove_duplicates(result, c("TF", "effect", "target"))
    if (nrow(result) > nrow(result_unique)) {
      warning(paste("Duplicates found in GRNI result for", file_name))
    }
    
    write.csv(result_unique, output_file, row.names = FALSE)
    
  }, error = function(e) {
    print(paste('Failure processing', bool_file, ':', e$message))
    fail_output <- file.path(output_folder, paste0('FAILGUROBI_', file_name, '.csv'))
    write.csv(bool_tfonly, fail_output, row.names = FALSE, quote = FALSE)
  })
}

print(paste("GRNI process completed. Results saved in:", output_folder))


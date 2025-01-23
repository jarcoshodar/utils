library(dplyr)
library(purrr)
library(cellxgenedp)


options(timeout = 36000)  #extremely high timeout so that no file is ever downloaded in half

#I created this with db <- db(), the cellxgenedp function of the library above
#then saveRDS in the folder this is launched from as db.rds. 
#this is to only download updates easily if we come back later by finding diff
#between this db and future ones
db <- readRDS("db.rds")

datasets <- datasets(db)
all_files <- files(db)

#we only want human data
human_datasets <- datasets %>%
  filter(map_lgl(organism, ~ any(map_chr(., "label") == "Homo sapiens")))

human_files <- all_files %>%
  filter(dataset_id %in% human_datasets$dataset_id)

#change if you want it elsewhere, in my case big data is read only on jobs
dest_dir <- "." 

if (!dir.exists(dest_dir)) {
  dir.create(dest_dir, recursive = TRUE)
}

#avoid confusing serial numbers
#can be matched back if you want some other filepart to be name
generate_filename <- function(file_info, extension) {
  dataset_id <- file_info$dataset_id
  title <- human_datasets %>%
    filter(dataset_id == file_info$dataset_id) %>%
    pull(title) %>%
    gsub("[^a-zA-Z0-9_]", "_", .)  # Replace non-alphanumeric characters with underscores
  
  paste0(title, "_", dataset_id, ".", extension)
}


human_files %>%
  rowwise() %>%
  mutate(download_status = {
    tryCatch({
      #as above
      extension <- tools::file_ext(url)
      new_filename <- generate_filename(cur_data(), extension)
      file_path <- file.path(dest_dir, new_filename)
      
      #skip .rds if .h5ad exists
      if (extension == "rds") {
        corresponding_h5ad <- sub("\\.rds$", ".h5ad", file_path)
        if (file.exists(corresponding_h5ad)) {
          return("Skipped: .rds skipped as .h5ad exists")
        }
      }
      
      #skip if file already exists
      if (file.exists(file_path)) {
        return("Skipped: File already exists")
      }
      
      #if none applies
      download.file(url, file_path, mode = "wb")
      "Success"
    }, error = function(e) {
      paste("Failed:", e$message)
    })
  })
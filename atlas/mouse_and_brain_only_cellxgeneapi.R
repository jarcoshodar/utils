library(dplyr)
library(tidyr)
library(purrr)
library(ontologyIndex)
library(cellxgenedp)
library(tools)

options(timeout = 36000)  #extremely high timeout so that no file is ever downloaded in half
# db <- db()
# saveRDS(db,"~/clean/atlas/db.rds")
#I created this with db <- db(), the cellxgenedp function of the library above
#then saveRDS in the folder this is launched from as db.rds. 
#this is to only download updates easily if we come back later by finding diff
#between this db and future ones
db <- readRDS("db.rds")

datasets <- datasets(db)
all_files <- files(db)

datasets <- datasets %>%
  filter(map_lgl(organism, ~ any(map_chr(., "label") == "Mus musculus")))

unnested <- datasets %>%
  select(dataset_id, tissue) %>%
  unnest_longer(tissue) %>%
  unnest_wider(tissue)

# head(unnested)
# 
# unique_labels <- unnested %>%
#   pull(label) %>%
#   unique() %>%
#   sort()
# 
# head(unique_labels, 20)

#approach to this manually was decided not to be viable so we use ontology

# wget http://purl.obolibrary.org/obo/uberon.obo 
#  uberon <- get_ontology("~/clean/atlas/uberon.obo", propagate_relationships = "is_a")
# fails in the full dataset to get that relationship. so we get:
# wget http://purl.obolibrary.org/obo/uberon/basic.obo
# 
uberon <- get_ontology(
  "basic.obo", 
  propagate_relationships = c("is_a", "part_of")
)

brain_descs <- get_descendants(uberon, "UBERON:0000955")

df_ont_based <- unnested %>%
  filter(ontology_term_id %in% brain_descs)

brain_dataset_ids <- unique(df_ont_based$dataset_id)

brain_files <- all_files %>%
  filter(dataset_id %in% brain_dataset_ids)

generate_filename <- function(file_info, extension, parent_datasets) {
  dataset_id <- file_info$dataset_id
  # find the dataset's 'title'
  title <- parent_datasets %>%
    filter(dataset_id == file_info$dataset_id) %>%
    pull(title) %>%
    gsub("[^a-zA-Z0-9_]", "_", .)  # Replace non-alphanumeric with underscores
  
  paste0(title, "_", dataset_id, ".", extension)
}

brain_datasets_info <- datasets %>% 
  select(dataset_id, title)

brain_files %>%
  rowwise() %>%
  mutate(download_status = {
    tryCatch({
      extension <- file_ext(url)
      new_filename <- generate_filename(cur_data(), extension, brain_datasets_info)
      file_path <- file.path(".", new_filename)
      
      # skip .rds if .h5ad exists
      if (extension == "rds") {
        corresponding_h5ad <- sub("\\.rds$", ".h5ad", file_path)
        if (file.exists(corresponding_h5ad)) {
          return("Skipped: .rds skipped as .h5ad exists")
        }
      }
      
      # skip if file already exists
      if (file.exists(file_path)) {
        return("Skipped: File already exists")
      }
      
      # if none applies, do the download
      download.file(url, file_path, mode = "wb")
      "Success"
    }, error = function(e) {
      paste("Failed:", e$message)
    })
  }) -> result_tibble


cat("\nDownload attempt summary:\n")
print(result_tibble %>% select(dataset_id, url, download_status))




#' Transfer cell type labels between Seurat objects
#' @param source_object Seurat object containing the cell type labels to transfer
#' @param target_object Seurat object to receive the labels
#' @param source_column Name of the column in source object containing cell types
#' @param target_column Name for the new column in target object
#' @return Modified target Seurat object with new cell type column
transfer_cell_labels <- function(source_object, 
                                 target_object, 
                                 source_column = "cell_type",
                                 target_column = "cell_type") {
  
  #get cell barcodes from both objects
  source_barcodes <- rownames(source_object@meta.data)
  target_barcodes <- rownames(target_object@meta.data)
  
  #create mapping from source object
  label_map <- source_object@meta.data[[source_column]]
  names(label_map) <- source_barcodes
  
  #create new column in target object
  target_object@meta.data[[target_column]] <- 
    #for each cell in target, either get its label from source or "no_data"
    ifelse(target_barcodes %in% source_barcodes,
           label_map[target_barcodes],
           "no_data")
  
  return(target_object)
}

best_seur <- transfer_cell_labels(
  source_object = old_seur,
  target_object = best_seur,
  source_column = "cell_type",  #column containing cell types in old_seur
  target_column = "cell_type"   #name for new column in best_seur
)

#preveiw
table(best_seur$cell_type, useNA = "ifany")
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

#gunction to integrate with specified reference (or NULL for no single dataset is used as referece)
integrate_seurat_pair <- function(file1, file2, reference = NULL, results_dir) {
  cat("Loading Seurat objects...\n")
  
  seurat_obj1 <- readRDS(file1)
  seurat_obj2 <- readRDS(file2)
  
  seurat_obj1$source_file <- basename(file1)
  seurat_obj2$source_file <- basename(file2)
  
  #ensure SCT assay
  if (!"SCT" %in% Assays(seurat_obj1) || !"SCT" %in% Assays(seurat_obj2)) {
    stop("Both objects must have SCT assay")
  }
  DefaultAssay(seurat_obj1) <- "SCT"
  DefaultAssay(seurat_obj2) <- "SCT"
  
  seurat_objects <- list(seurat_obj1, seurat_obj2)
  

  features <- SelectIntegrationFeatures(object.list = seurat_objects, nfeatures = 6000)
  seurat_objects <- PrepSCTIntegration(object.list = seurat_objects, 
                                       anchor.features = features)
  

  ref_str <- if(is.null(reference)) "no_ref" else paste0("ref_", reference)
  cat("Finding anchors with", ref_str, "...\n")
  
  anchors <- FindIntegrationAnchors(object.list = seurat_objects,
                                    normalization.method = "SCT",
                                    anchor.features = features,
                                    reference = reference)
  

  integrated <- IntegrateData(anchorset = anchors,
                              normalization.method = "SCT")
  

  output_file <- sprintf("integrated_%s.rds", ref_str)
  saveRDS(integrated, file.path(results_dir, output_file))
  cat("Integration complete. Result saved as", output_file, "\n")
  
  return(integrated)
}


if (!dir.exists("results")) dir.create("results")

cat("\nTrying with first object as reference...\n")
integrate_seurat_pair("SeuratObject_p2.rds", "SeuratObject_p3.rds", 
                      reference = 1, results_dir = "results")

cat("\nTrying with second object as reference...\n")
integrate_seurat_pair("SeuratObject_p2.rds", "SeuratObject_p3.rds", 
                      reference = 2, results_dir = "results")

cat("\nTrying without reference...\n")
integrate_seurat_pair("SeuratObject_p2.rds", "SeuratObject_p3.rds", 
                      reference = NULL, results_dir = "results")

cat("All integration attempts completed.\n")


#!/usr/bin/env Rscript
numPCs <- 50
umap_neighbors <- 20

suppressPackageStartupMessages({
  library(Seurat)           
  library(monocle3)         
  library(SeuratWrappers)  
  library(ggplot2)
  library(fpc)             
  library(mclust)           # for ARI
  library(entropy)          
  library(lisi)             
})
#nokbet due to computational demandiness, if ever in separate script
#adjust these paths as needed:
input_dir  <- "~/clean/sybille/metainvivo_out_s4_all_5_more_genes/"
output_dir <- "~/clean/sybille/rate_more_genes/"

#expand tilde (~) if present
input_dir  <- path.expand(input_dir)
output_dir <- path.expand(output_dir)

if (!dir.exists(input_dir)) {
  stop("Input directory does not exist: ", input_dir)
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive=TRUE)
  message("Created output directory: ", output_dir)
}

resolutions <- c(rbind(10^seq(-7,-5), 5 * 10^seq(-7,-5)))
# Example expansions: 1e-07, 5e-07, 1e-06, 5e-06, 1e-05, 5e-05

#assume your dataset column is named "source_file"
dataset_col <- "source_file"


script_start <- Sys.time()
cat("Script started at:", format(script_start), "\n")

computeCH <- function(emb_pca, partition_vec) {
  fpc::calinhara(emb_pca, as.integer(partition_vec))
}

computeARI <- function(partition_vec, batch_vec) {
  cluster_labels <- as.numeric(as.factor(partition_vec))
  dataset_labels <- as.numeric(as.factor(batch_vec))
  mclust::adjustedRandIndex(cluster_labels, dataset_labels)
}

computeEntropy <- function(partition_vec, batch_vec) {
  cluster_levels <- unique(partition_vec)
  ent_vals <- sapply(cluster_levels, function(cl) {
    cells_idx <- which(partition_vec == cl)
    ds_counts <- table(batch_vec[cells_idx])
    ds_probs  <- ds_counts / sum(ds_counts)
    -sum(ds_probs * log(ds_probs))  # Shannon entropy
  })
  data.frame(Cluster=cluster_levels, Entropy=ent_vals)
}

computeClusterDatasetStats <- function(partition_vec, batch_vec) {
  ct <- table(partition_vec, batch_vec)
  prop <- prop.table(ct, margin = 1)
  list(ct_table = ct, ct_prop = prop)
}

computeLISI <- function(seur, partition_vec, batch_col="source_file", perplexity=15) {
  if (!"umap" %in% names(seur@reductions)) {
    warning("No UMAP found in the Seurat object for LISI. Returning NA.")
    return(NA)
  }
  
  umap_emb <- Embeddings(seur, "umap")
  meta_df  <- data.frame(
    batch = seur@meta.data[[batch_col]],
    row.names = colnames(seur)
  )
  
  #check rownames alignment
  cat("UMAP rownames length:", length(rownames(umap_emb)), "\n")
  cat("meta_df rownames length:", length(rownames(meta_df)), "\n")
  
  all_same <- all(rownames(umap_emb) == rownames(meta_df))
  cat("All rownames identical? ", all_same, "\n")
  
  if (!all_same) {
    common_cells <- intersect(rownames(umap_emb), rownames(meta_df))
    cat("Number of common cells:", length(common_cells), "\n")
    
    if (length(common_cells) == 0) {
      warning("No common cells found between embeddings and metadata. Skipping LISI.")
      return(NA)
    }
    
    umap_emb <- umap_emb[common_cells, , drop = FALSE]
    meta_df  <- meta_df[common_cells, , drop = FALSE]
    cat("Subsetting done. New rownames length:", length(rownames(umap_emb)), "\n")
  }
  
  na_emb <- sum(is.na(umap_emb))
  if (na_emb > 0) {
    cat("Number of NA values in UMAP embeddings:", na_emb, "\n")
    warning("UMAP embeddings contain NA values. Consider removing these cells.")
    return(NA)
  }
  #we treat experimental source as batch
  num_batches <- length(unique(meta_df$batch))
  cat("Number of unique batches:", num_batches, "\n")
  
  batch_counts <- table(meta_df$batch)
  print(head(batch_counts))
  
  if (any(batch_counts < (3 * perplexity + 1))) {
    warning("Some batches have fewer cells than required for LISI. Consider adjusting perplexity or excluding small batches.")
  }
  
  lisi_res <- tryCatch({
    lisi::compute_lisi(
      X = umap_emb,
      meta_data = meta_df,
      label_colnames = "batch",
      perplexity = perplexity
    )
  }, error = function(e) {
    warning("LISI computation failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(lisi_res)) {
    return(NA)
  }
  
  
  na_lisi <- sum(is.na(lisi_res$batch))
  if (na_lisi > 0) {
    cat("Number of NA LISI scores:", na_lisi, "\n")
  }
  
  #calculate mean LISI score for the 'batch' label
  mean_lisi <- mean(lisi_res$batch, na.rm=TRUE)
  cat("Mean LISI score for 'batch': ", mean_lisi, "\n")
  
  return(mean_lisi)
}


#supposed to be different integrations to compare
all_rds_files <- list.files(
  path        = input_dir,
  pattern     = "^integrated_object_ref_.*\\.rds$",
  full.names  = TRUE
)
if (length(all_rds_files) == 0) {
  stop("No matching RDS files found in input directory.")
}

for (rds in all_rds_files) {
  file_start <- Sys.time()
  cat("\n>>> Processing:", rds, "started at", format(file_start), "\n")
  
  base_name <- tools::file_path_sans_ext(basename(rds))

  cat("[", format(Sys.time()), "] Loading Seurat object...\n")
  seur <- readRDS(rds)
  cat("[", format(Sys.time()), "] Seurat object loaded with",
      ncol(seur), "cells and", nrow(seur), "features.\n")
  
  # --(A) Make sure we have PCA--
  if (!"pca" %in% names(seur@reductions)) {
    seur <- RunPCA(seur, npcs = 100, verbose = TRUE)
  }
  cat("[", format(Sys.time()), "] Extracting PCA embeddings...\n")
  emb_pca <- Embeddings(seur, "pca")  # Use all PCs by default
  
  # --(B) Run Seurat's FindNeighbors for consistency--
  cat("[", format(Sys.time()), "] Running Seurat::FindNeighbors...\n")
  seur <- FindNeighbors(seur, dims = 1:numPCs, verbose=FALSE)
  cat("[", format(Sys.time()), "] Done with FindNeighbors.\n")
  
  # --(C) We won't do Seurat's FindClusters here because we are relying on Monocle,
  #        but if you prefer a default Seurat cluster to check for NAs, you could do it.
  
  # --(D) Run UMAP if missing--
  if (!"umap" %in% names(seur@reductions)) {
    cat("[", format(Sys.time()), "] Running UMAP (Seurat)...\n")
    seur <- RunUMAP(seur, dims = 1:numPCs, n.neighbors = umap_neighbors, verbose=FALSE)
    cat("[", format(Sys.time()), "] UMAP done.\n")
  } else {
    cat("[", format(Sys.time()), "] UMAP already found in object. Not re-running.\n")
  }
  
  #build a base CellDataSet manually (necessary because of version differences)
  cat("[", format(Sys.time()), "] Preparing expression matrix and metadata...\n")
  expression_matrix <- GetAssayData(seur, assay = "RNA", slot = "data")
  cell_metadata     <- seur@meta.data
  gene_metadata     <- data.frame(
    gene_short_name = rownames(expression_matrix),
    row.names       = rownames(expression_matrix)
  )
  cat("[", format(Sys.time()), "] Expression matrix and metadata prepared.\n")
  
  cat("[", format(Sys.time()), "] Creating base CellDataSet object...\n")
  base_cds <- new_cell_data_set(
    expression_data = expression_matrix,
    cell_metadata   = cell_metadata,
    gene_metadata   = gene_metadata
  )
  cat("[", format(Sys.time()), "] Base CellDataSet object created.\n")
  
  #add UMAP reductions to the CDS if they exist
  if (is.null(reducedDims(base_cds))) {
    reducedDims(base_cds) <- S4Vectors::SimpleList()
  }
  if ("umap" %in% names(seur@reductions)) {
    umap_embeddings <- Embeddings(seur, "umap")
    reducedDims(base_cds)$UMAP <- umap_embeddings
  } else {
    warning("No UMAP found in Seurat object. cluster_cells(UMAP) may fail in Monocle.")
  }
  

  results_for_this_object <- list()
  dataset_values <- seur@meta.data[[dataset_col]]
  
  #track CH index for each resolution to pick a "winning" resolution
  ch_vals <- numeric(length(resolutions))
  names(ch_vals) <- as.character(resolutions)
  

  for (i in seq_along(resolutions)) {
    res <- resolutions[i]
    res_start <- Sys.time()
    cat("\n[", format(res_start), "] Starting Monocle clustering at resolution =", res, "\n")
    

    cds_local <- base_cds
    

    cat("[", format(Sys.time()), "] cluster_cells (Leiden) in Monocle3...\n")
    cds_local <- tryCatch({
      cluster_cells(
        cds_local,
        k                = 10,
        reduction_method = "UMAP",
        cluster_method   = "leiden",
        resolution       = res,
        num_iter         = 5
      )
    }, error = function(e) {
      warning("cluster_cells failed for resolution=", res, ": ", e$message)
      cat("[", format(Sys.time()), "] Clustering failed for resolution =", res, "\n")
      return(NULL)
    })
    
    if (is.null(cds_local)) {
      cat("[", format(Sys.time()), "] Skipping resolution =", res, "due to error.\n")
      next
    }
    cat("[", format(Sys.time()), "] Clustering completed for resolution =", res, "\n")
    

    part_vec <- monocle3::clusters(cds_local)
    num_clusters <- length(unique(part_vec))
    cat("[", format(Sys.time()), "] Number of clusters formed:", num_clusters, "\n")
    

    na_count <- sum(is.na(part_vec))
    cat("      -> # of NA cluster assignments for this resolution:", na_count,
        "out of", length(part_vec), "cells.\n")
    
    if (num_clusters > 50) {
      cat("[", format(Sys.time()), "] SKIPPING resolution =", res, 
          "because it yielded", num_clusters, "clusters ( > 50 )\n")
      rm(cds_local)
      gc()
      next
    }
    
    cat("[", format(Sys.time()), "] Computing Calinski–Harabasz index...\n")
    ch_val <- computeCH(emb_pca, part_vec)
    ch_vals[as.character(res)] <- ch_val
    
    cat("[", format(Sys.time()), "] Computing Adjusted Rand Index (ARI)...\n")
    ari_val <- computeARI(part_vec, dataset_values)
    
    cat("[", format(Sys.time()), "] Computing cluster-dataset statistics...\n")
    cd_stats <- computeClusterDatasetStats(part_vec, dataset_values)
    
    cat("[", format(Sys.time()), "] Computing cluster-level entropy...\n")
    ent_df  <- computeEntropy(part_vec, dataset_values)
    avg_ent <- mean(ent_df$Entropy)
    
    cat("[", format(Sys.time()), "] Computing LISI score...\n")
    lisi_val <- NA
    tryCatch({
      lisi_val <- computeLISI(seur, part_vec, batch_col = dataset_col, perplexity = 15)
    }, error = function(e) {
      warning("LISI error at resolution=", res, ": ", e$message)
      cat("[", format(Sys.time()), "] LISI computation failed for resolution =", res, "\n")
    })
    
    outprefix <- file.path(output_dir, paste0(base_name, "_res", res))
    write.csv(as.data.frame(cd_stats$ct_table),
              file = paste0(outprefix, "_cluster_dataset_table.csv"),
              row.names=TRUE)
    write.csv(round(cd_stats$ct_prop, 2),
              file = paste0(outprefix, "_cluster_dataset_prop.csv"),
              row.names=TRUE)
    write.csv(ent_df,
              file=paste0(outprefix, "_entropy.csv"),
              row.names=FALSE)
    
    #summarize row for this resolution
    row_df <- data.frame(
      Resolution        = res,
      NumClusters       = num_clusters,
      ARI               = ari_val,
      MeanEntropy       = avg_ent,
      MeanLISI          = lisi_val,
      CalinskiHarabasz  = ch_val
    )
    results_for_this_object[[as.character(res)]] <- row_df
    

    rm(cds_local)
    gc()
    
    res_end <- Sys.time()
    cat("[", format(res_end), "] Completed resolution =", res, "in",
        round(difftime(res_end, res_start, units="secs"), 2), "seconds.\n")
  }  
  
  #combine the rows for all resolutions
  if (length(results_for_this_object) > 0) {
    cat("[", format(Sys.time()), "] Aggregating results for all resolutions...\n")
    combined_df <- do.call(rbind, results_for_this_object):
    out_agg <- file.path(output_dir, paste0(base_name, "_allRes_metrics.csv"))
    write.csv(combined_df, file=out_agg, row.names=FALSE)
    cat("[", format(Sys.time()), "] Saved aggregator CSV:", out_agg, "\n")
    
    # (E) Pick the "winning" resolution by highest CH index
    best_res <- names(which.max(ch_vals))
    cat("[", format(Sys.time()), "] Best resolution by CH index is:", best_res, "\n")
    
    #re-run cluster_cells at that best resolution
    cat("[", format(Sys.time()), "] Re-clustering at best resolution to assign final labels...\n")
    cds_final <- base_cds
    cds_final <- cluster_cells(
      cds_final,
      k                = 10,
      reduction_method = "UMAP",
      cluster_method   = "leiden",
      resolution       = as.numeric(best_res),
      num_iter         = 5
    )
    final_part_vec <- clusters(cds_final)
    cat(" -> #clusters:", length(unique(final_part_vec)), "\n")
    cat(" -> #NAs:", sum(is.na(final_part_vec)), "\n")
    
    #assign these final Monocle clusters to the Seurat object:
    seur$seurat_clusters <- final_part_vec
    Idents(seur)         <- final_part_vec
    
    # (Optionally) you could save the final updated Seurat object
    #change below if you want to have same filename (overwrite)
    #saveRDS(seur, file.path(output_dir, paste0(base_name, "_final_withBestRes.rds")))
  } else {
    cat("[", format(Sys.time()), "] No resolutions processed for", base_name, " (list was empty?)\n")
  }
  
  #time for this file
  file_end <- Sys.time()
  cat(">>> Finished:", rds, "at", format(file_end), "\n")
  cat("Time elapsed for this file:",
      round(difftime(file_end, file_start, units="secs"), 2),
      "seconds\n")
}
#timeforall
script_end <- Sys.time()
cat("\nAll files processed.\n")
cat("Total script runtime:",
    round(difftime(script_end, script_start, units="secs"), 2),
    "seconds\n")
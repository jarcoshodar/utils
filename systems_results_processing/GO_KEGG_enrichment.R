#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(WebGestaltR)
  library(tools)  
})

option_list <- list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="Path to a single CSV file with DE results. If provided, overrides --folder and --root/--study mode."),
  make_option(c("--folder"), type="character", default=NULL,
              help="Path to a folder containing CSV files to process."),
  make_option(c("--root"), type="character", default=NULL,
              help="Root path containing <study>/results3/Integrated/<organ> subfolders with CSV files. Used if --file and --folder are not provided."),
  make_option(c("--study"), type="character", default=NULL,
              help="Study name under --root. Combined with --root to find CSV files in 'results3/Integrated/*'."),
  make_option(c("--output_dir"), type="character", default=NULL,
              help="Output directory where enrichment results will be saved."),
  make_option(c("--organism"), type="character", default="mmusculus",
              help="Organism for WebGestaltR (e.g., 'mmusculus' or 'rnorvegicus'). [default: %default]"),
  make_option(c("--pval_threshold"), type="double", default=0.05,
              help="P-value (or adjusted p-value) threshold to consider genes significant. [default: %default]"),
  make_option(c("--padj_column"), type="character", default="p_val",
              help="Column name in the CSV for the p-value (or adjusted p-value) to use for significance. [default: %default]"),
  make_option(c("--logfc_column"), type="character", default="avg_log2FC",
              help="Column name in the CSV for log2 fold-change. [default: %default]"),
  make_option(c("--run_go"), type="logical", default=TRUE,
              help="Whether to run GO enrichment (separately for up- and down-regulated genes). [default: %default]"),
  make_option(c("--run_kegg"), type="logical", default=TRUE,
              help="Whether to run KEGG enrichment (using combined significant genes). [default: %default]"),
  make_option(c("--background_file"), type="character", default=NULL,
              help="Path to a text file listing background genes (one per line)."),
  make_option(c("--seurat_file"), type="character", default=NULL,
              help="(Optional) Path to a Seurat .rds file from which to extract all rownames as background. If provided, overrides --background_file.")
)

parser <- OptionParser(option_list = option_list,
                       usage = "%prog [options]",
                       description = "Perform WebGestaltR enrichment analysis on CSV-based DE results.
By default, the script will search for a single CSV (--file), or a folder of CSVs (--folder),
or use a root/study approach (--root and --study).
For each CSV, significant genes are determined using --pval_threshold (from the column specified by --padj_column)
and split into up-/down-regulated sets for GO enrichment (if --run_go is TRUE),
while the combined set is used for KEGG enrichment (if --run_kegg is TRUE).
Results are saved in subfolders (pGO, nGO, KEGG) under the directory given by --output_dir."
)

opt <- parse_args(parser)


csv_files <- character()

if (!is.null(opt$file)) {
  #Single-file mode
  csv_files <- c(opt$file)
} else if (!is.null(opt$folder)) {
  if (!dir.exists(opt$folder)) {
    stop("Folder does not exist: ", opt$folder)
  }
  csv_files <- list.files(opt$folder, pattern="\\.csv$", full.names=TRUE)
} else if (!is.null(opt$root) && !is.null(opt$study)) {
  base_path <- file.path(opt$root, opt$study, "results3", "Integrated")
  if (!dir.exists(base_path)) {
    stop("Path not found: ", base_path)
  }
  subfolders <- list.dirs(base_path, recursive=FALSE)
  for (fld in subfolders) {
    these_csvs <- list.files(fld, pattern="\\.csv$", full.names=TRUE)
    csv_files <- c(csv_files, these_csvs)
  }
} else {
  cat("\nNo --file, no --folder, and no valid --root/--study given. Nothing to do.\n\n")
  print_help(parser)
  quit(save="no", status=1)
}

if (length(csv_files) == 0) {
  cat("No CSV files found to process. Exiting.\n")
  quit(save="no", status=0)
}

do_webgestalt <- function(organism, gene_list, background_genes, database, label, fdrThr=0.01) {
  if (length(gene_list) == 0) {
    cat("No genes in interestGene, skipping.\n")
    return(data.frame())
  }
  
  #Remove any NA or blank values
  gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]
  background_genes <- background_genes[!is.na(background_genes) & background_genes != ""]
  
  if (length(gene_list) == 0) {
    cat("All genes were NA or blank, skipping.\n")
    return(data.frame())
  }
  
  cat("WebGestaltR:", label, " => #interestGenes =", length(gene_list),
      " #background =", length(background_genes), "\n")
  
  res_df <- data.frame()
  tryCatch({
    res_df <- WebGestaltR(
      organism         = organism,
      enrichMethod     = "ORA",
      enrichDatabase   = database,
      interestGene     = gene_list,
      interestGeneType = "genesymbol",
      referenceGene    = background_genes,
      referenceGeneType= "genesymbol",
      minNum           = 10,
      maxNum           = 500,
      sigMethod        = "fdr",
      fdrMethod        = "BH",
      fdrThr           = fdrThr,
      topThr           = 10,
      isOutput         = FALSE,
      collapseMethod   = "mean"
    )
    cat("Success:", label, " => rows in result:", nrow(res_df), "\n")
  }, error = function(e) {
    cat("Error in WebGestaltR:", conditionMessage(e), "\n")
  })
  
  if (is.null(res_df) || nrow(res_df) == 0) {
    res_df <- data.frame()
  }
  res_df
}


process_one_csv <- function(csv_file, opt) {
  cat("\nProcessing CSV:", csv_file, "\n")
  
  df <- read.csv(csv_file, header=TRUE, stringsAsFactors=FALSE)
  
  #determine background genes:
  background_genes <- character()
  if (!is.null(opt$seurat_file)) {
    if (!file.exists(opt$seurat_file)) {
      stop("Seurat file not found:", opt$seurat_file)
    }
    cat("Reading Seurat object for background from:", opt$seurat_file, "\n")
    sobj <- readRDS(opt$seurat_file)
    background_genes <- rownames(sobj)
  } else if (!is.null(opt$background_file)) {
    if (!file.exists(opt$background_file)) {
      stop("Background file not found:", opt$background_file)
    }
    cat("Reading background from text file:", opt$background_file, "\n")
    background_genes <- readLines(opt$background_file)
  } else {
    stop("You must provide either --seurat_file or --background_file for background genes.")
  }
  

  if (!(opt$padj_column %in% colnames(df))) {
    stop("Specified p-value column (", opt$padj_column, ") not found in CSV:", csv_file)
  }
  if (!(opt$logfc_column %in% colnames(df))) {
    stop("Specified logFC column (", opt$logfc_column, ") not found in CSV:", csv_file)
  }
  
  sig_rows <- df[df[[opt$padj_column]] < opt$pval_threshold, ]
  if (nrow(sig_rows) == 0) {
    cat("No significant genes below threshold", opt$pval_threshold, " => no up/down sets.\n")
  }
  
  up_genes <- sig_rows[sig_rows[[opt$logfc_column]] > 0, "Gene"]
  down_genes <- sig_rows[sig_rows[[opt$logfc_column]] < 0, "Gene"]
  all_sig_genes <- sig_rows[["Gene"]]
  
  out_dir <- if (!is.null(opt$output_dir)) {
    opt$output_dir
  } else {
    dirname(csv_file)
  }
  

  save_enrichment <- function(result_df, label) {
    if (nrow(result_df) == 0) {
      out_file <- file.path(out_dir, paste0("FAILURE_", label, ".csv"))
    } else {
      out_file <- file.path(out_dir, paste0(label, ".csv"))
    }
    write.csv(result_df, out_file, row.names=FALSE)
    cat("Wrote:", out_file, "\n")
  }
  
  #GO enrichment (up and down separately)
  if (opt$run_go) {
    pGO_folder <- file.path(out_dir, "pGO")
    nGO_folder <- file.path(out_dir, "nGO")
    if (!dir.exists(pGO_folder)) dir.create(pGO_folder, recursive=TRUE)
    if (!dir.exists(nGO_folder)) dir.create(nGO_folder, recursive=TRUE)
    
    if (length(up_genes) > 0) {
      label_up <- paste0(file_path_sans_ext(basename(csv_file)), "_Positive_GO_Enrichment")
      res_up <- do_webgestalt(
        organism = opt$organism, 
        gene_list = up_genes,
        background_genes = background_genes,
        database = "geneontology_Biological_Process",
        label = label_up,
        fdrThr = 0.01
      )
      save_enrichment(res_up, label_up)
    }
    
    if (length(down_genes) > 0) {
      label_down <- paste0(file_path_sans_ext(basename(csv_file)), "_Negative_GO_Enrichment")
      res_down <- do_webgestalt(
        organism = opt$organism,
        gene_list = down_genes,
        background_genes = background_genes,
        database = "geneontology_Biological_Process",
        label = label_down,
        fdrThr = 0.01
      )
      save_enrichment(res_down, label_down)
    }
  }
  
  #KEGG enrichment (combined up+down)
  if (opt$run_kegg) {
    kegg_folder <- file.path(out_dir, "KEGG")
    if (!dir.exists(kegg_folder)) dir.create(kegg_folder, recursive=TRUE)
    
    label_kegg <- paste0(file_path_sans_ext(basename(csv_file)), "_KEGG_Enrichment")
    res_kegg <- do_webgestalt(
      organism = opt$organism,
      gene_list = all_sig_genes,
      background_genes = background_genes,
      database = "pathway_KEGG",
      label = label_kegg,
      fdrThr = 0.01
    )
    save_enrichment(res_kegg, label_kegg)
  }
}


for (f in csv_files) {
  if (!file.exists(f)) {
    cat("File does not exist, skipping:", f, "\n")
    next
  }
  process_one_csv(f, opt)
}

cat("\nDone. Processed", length(csv_files), "CSV files.\n")

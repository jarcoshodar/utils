#!/usr/bin/env Rscript
# deldge_opt.R
#
#  1) Reads single or multiple Seurat .rds files
#  2) Expects a "cell type" column and a "condition" column in meta.data
#  3) Creates a combined 'cellcon' = paste0(celltype, condition)
#  4) For each cell type, runs pseudobulk dge between e.g. "McrCR" vs "McrO" (if cond1=CR, cond2=O)
#  5) Writes one .csv per cell type, plus one background file per .rds
## NOT IMPLEMENTED YET: use of 'real replicates, if given'
## NO IMPLEMENTED YET: able to adjust for confounders in a deseq ~ + fashion (need to modify delegate library)
suppressPackageStartupMessages({
  library(optparse)   
  library(Seurat)
  library(DElegate)    #pseudobulk. see https://www.biorxiv.org/content/10.1101/2023.03.28.534443v2
  library(stringr)
})
# Expected Input Structure:
# /path/to/root/
# └── StudyName/
#     └── results/
#         ├── Integrated_Organ1.rds       
#         ├── Integrated_Organ2.rds
#         └── Integrated_Organ3.rds
#
# Output Structure:
# /path/to/root/
# └── StudyName/
#     └── results/
#         ├── Integrated/                
#         │   ├── Organ1/                
#         │   │   ├── CellType1Cond1_vs_CellType1Cond2.csv
#         │   │   ├── CellType2Cond1_vs_CellType2Cond2.csv
#         │   │   └── background_Cond1_vs_Cond2.txt
#         │   ├── Organ2/
#         │   │   ├── CellType1Cond1_vs_CellType1Cond2.csv
#         │   │   └── background_Cond1_vs_Cond2.txt
#         │   └── Organ3/
#         │       ├── CellType1Cond1_vs_CellType1Cond2.csv
#         │       └── background_Cond1_vs_Cond2.txt
#         ├── Integrated_Organ1.rds       
#         ├── Integrated_Organ2.rds
#         └── Integrated_Organ3.rds
option_list <- list(
  make_option(
    c("--file"), type="character", default=NULL,
    help="Path to a single Seurat .rds file. Overrides --root & --study if given."
  ),
  make_option(
    c("--root"), type="character", default=NULL,
    help="Root path containing the study folder with results/Integrated_*.rds files."
  ),
  make_option(
    c("--study"), type="character", default=NULL,
    help="Study name under 'root'. E.g. 'Parabiosis'."
  ),
  make_option(
    c("--organs"), type="character", default=NULL,
    help="Comma-separated organ names, e.g. 'Aorta,Cerebellum'. Looks for corresponding Integrated_*.rds files in (root/study)/results/. If omitted, all Integrated_*.rds files in results/ are processed."
  ),

  make_option(
    c("--celltype_col"), type="character", default="cell_type",
    help="Name of meta.data column with cell type labels (default: '%default')."
  ),
  make_option(
    c("--condition_col"), type="character", default="condition",
    help="Name of meta.data column with condition labels (default: '%default')."
  ),
  make_option(
    c("--condition1"), type="character", default=NULL,
    help="First value of condition_col (ident.1) to compare. REQUIRED."
  ),
  make_option(
    c("--condition2"), type="character", default=NULL,
    help="Second value of condition_col label (ident.2) to compare. Results are thus 'condition1 over condition2' (positive log fold change means more expressed in condition1) REQUIRED."
  ),

  make_option(
    c("--filtered"), type="logical", default=TRUE,
    help="Whether to filter results to keep only genes with |log2FC|>=0.25 and expressed in >10%% of cells in at least one condition. [default: %default]."
  ),
  make_option(
    c("--background"), type="logical", default=TRUE,
    help="Whether to save one background .txt per .rds (listing all measured genes). [default: %default]."
  )
)
parser <- OptionParser(
  usage = "%prog [options]",
  description = "Perform cell-type + condition (cellcon) DE. 
For each .rds, creates a combined 'cellcon' = paste0(celltype, condition) 
and compares e.g. 'McrCR' vs 'McrO' for each cell type. 
Writes one CSV per cell type, plus a single background if requested. Background is just genes measured, useful for downstream e.g. WebGestaltR"
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$condition1) || is.null(opt$condition2)) {
  cat("\nERROR: You must provide --condition1 and --condition2.\n\n")
  print_help(parser)
  quit(save="no", status=1)
}
file_list <- NULL
if (!is.null(opt$file)) {
  #single-file mode
  file_list <- list(opt$file)
} else {
  #multi-file mode => need --root & --study
  if (is.null(opt$root) || is.null(opt$study)) {
    cat("\nERROR: Either provide --file or use --root and --study.\n\n")
    print_help(parser)
    quit(save="no", status=1)
  }
  #if organs provided, build path for each organ
  if (!is.null(opt$organs)) {
    orgs <- unlist(strsplit(opt$organs, ","))
    file_list <- file.path(opt$root, opt$study, "results", paste0("Integrated_", orgs, ".rds"))  } else {
    #else gather all rds in results/ that match pattern
    rds_files <- list.files(
      path       = file.path(opt$root, opt$study, "results"),
      pattern    = "^Integrated_.*\\.rds$",      
      full.names = TRUE
    )
    file_list <- rds_files
  }
}
if (length(file_list) == 0) {
  cat("\nNo RDS files found to process. Exiting.\n")
  quit(save="no", status=0)
}

run_dge_cellcon <- function(rds_path,
                            celltype_col,
                            condition_col,
                            cond1,
                            cond2,
                            do_filter,
                            do_background) {

  cat("\nProcessing:", rds_path, "\n")
  #parse organ from filename "Integrated_<organ>.rds"
  base_name <- basename(rds_path)
  organ <- sub("^Integrated_(.*)\\.rds$", "\\1", base_name)
  #define output dir => e.g. /.../study/results/Integrated/<organ>/
  top_dir <- dirname(dirname(rds_path))  # e.g. /.../study
  out_dir <- file.path(top_dir, "results", "Integrated", organ)
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

  #read object
  seur <- readRDS(rds_path)
  DefaultAssay(seur) <- "RNA"  #or "SCT" but deseq expects raw counts

  #check that celltype_col & condition_col exist
  meta_df <- seur@meta.data
  if (!(celltype_col %in% colnames(meta_df))) {
    stop("Column '", celltype_col, "' not found in Seurat metadata.")
  }
  if (!(condition_col %in% colnames(meta_df))) {
    stop("Column '", condition_col, "' not found in Seurat metadata.")
  }

  #subset the object to only cond1 & cond2
  keep_idx <- which(meta_df[[condition_col]] %in% c(cond1, cond2))
  if (length(keep_idx) < 2) {
    stop("No cells found for conditions (", cond1, " or ", cond2, ") in '", condition_col, "' column.")
  }
  seur <- subset(seur, cells=rownames(meta_df)[keep_idx])

  #create new 'cellcon' col => paste0(celltype_val, condition_val)
  #e.g. if celltype_col="cell_type" => "Mcr"; condition_col="condition" => "CR"
  #=> "McrCR"
  seur@meta.data$cellcon <- paste0(seur@meta.data[[celltype_col]],
                                   seur@meta.data[[condition_col]])

  #set Idents(seur) to that combined col
  Idents(seur) <- seur@meta.data$cellcon

  #get unique cell types from the cell_type column
  #for each cell type t, we expect t+cond1 and t+cond2 
  #e.g. "McrCR" vs "McrO"
  all_celltypes <- unique(seur@meta.data[[celltype_col]])

  #We'll keep track of how many comparisons we actually performed
  comp_count <- 0

  for (ct in all_celltypes) {
    group1 <- paste0(ct, cond1)
    group2 <- paste0(ct, cond2)

    #check if both groups exist
    these_levels <- levels(seur)  # all Idents
    if ((group1 %in% these_levels) && (group2 %in% these_levels)) {

      cat("Running DE for cell type:", ct, "=>", group1, "vs", group2, "\n")
  
      res <- findDE(seur, compare=c(group1, group2), method="deseq")

      #keep columns consisent with deseq/find markers
      colnames(res)[colnames(res)=="pvalue"]  <- "p_val"
      colnames(res)[colnames(res)=="log_fc"]  <- "avg_log2FC"
      colnames(res)[colnames(res)=="rate1"]   <- "pct.1"
      colnames(res)[colnames(res)=="rate2"]   <- "pct.2"
      colnames(res)[colnames(res)=="padj"]    <- "p_val_adj"
      colnames(res)[colnames(res)=="feature"] <- "Gene"

      
      if (do_filter) {
        res <- res[ abs(res$avg_log2FC)>=0.25 & (res$pct.1>=0.1 | res$pct.2>=0.1), ]
      }

    
      out_file <- file.path(out_dir, paste0(group1, "_vs_", group2, ".csv"))    
      cat("Saving DGE results to:", out_file, "\n")
      write.csv(res, out_file, row.names=FALSE, quote=FALSE)

      comp_count <- comp_count + 1
    }
  } 

  cat("Total comparisons done for", rds_path, ":", comp_count, "\n")

  if (do_background) {
    bg_file <- file.path(out_dir, paste0("background_", cond1, "_vs_", cond2, ".txt"))
    cat("Saving background genes to:", bg_file, "\n")
    writeLines(rownames(seur), con=bg_file)
  }
}

for (f in file_list) {
  run_dge_cellcon(
    rds_path      = f,
    celltype_col  = opt$celltype_col,
    condition_col = opt$condition_col,
    cond1         = opt$condition1,
    cond2         = opt$condition2,
    do_filter     = opt$filtered,
    do_background = opt$background
  )
}
cat("\nAll done.\n")
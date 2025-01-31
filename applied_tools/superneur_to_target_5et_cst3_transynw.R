#!/usr/bin/env Rscript

setwd('~/clean/transsynw/')
## Usage: Rscript code/SynergisticCore_noExpCutoff_sizeNormalization.R sample_data_Gokce2016.csv ::Oligo::Astro:: . sample_data_Gokce2016_cluster.csv hESC_E-MTAB-6819_H9_naive_smartseq2.Robj Mouse sample_data::sample_cluster
simulateArgs <- function(cmdLine) {
  # Simulate splitting the command line as if entered in a shell
  # Split on spaces to mimic shell argument parsing
  allParts <- strsplit(cmdLine, " ")[[1]]
  
  # Ignore the first two elements simulating the Rscript call and script name
  args <- allParts[-(1:2)]
  
  # Return the args vector
  return(args)
}

# Example command line
#cmdLine <- "Rscript code/SynergisticCore_noExpCutoff_sizeNormalization.R /home/jarcos/clean/sisifo/for_transsynw/cst3/all_combined_gene_expression_matrix.tsv ::L5_ET:: cst3_ast_to_super_neuron/ /home/jarcos/clean/sisifo/for_transsynw/cst3/nncell_annotation.tsv mAstro_GSE114000_smartseq2.Robj Mouse cst3::el5"

#cmdLine2 <- "Rscript code/SynergisticCore_noExpCutoff_sizeNormalization.R /storage/sisifo/for_transsynw/cst3/all_combined_gene_expression_matrix.tsv ::L5_ET:: cst3ownastafternn/ /storage/sisifo/for_transsynw/cst3/nncell_annotation.tsv /storage/sisifo/for_transsynw/cst3/astro_gene_expression_matrix.tsv Mouse cst3::el5"
#cmdLine22 <- "Rscript code/SynergisticCore_noExpCutoff_sizeNormalization.R /storage/sisifo/for_transsynw/cst3/all_combined_gene_expression_matrix.tsv ::L5_ET:: l5means_neuron_no_noise /storage/sisifo/for_transsynw/cst3/nncell_annotation.tsv /storage/sisifo/for_transsynw/cst3/astro_gene_expression_matrix.tsv Mouse cst3::el5"
cmdLine22 <- "Rscript code/SynergisticCore_noExpCutoff_sizeNormalization.R /home/jarcos/clean/sisifo/for_transsynw/cst3/all_combined_gene_expression_matrix.tsv ::L5_ET:: from_other_N_to5ET /home/jarcos/clean/sisifo/for_transsynw/cst3/nncell_annotation.tsv /home/jarcos/clean/sisifo/for_transsynw/cst3/otherN_gene_expression_matrix.tsv Mouse cst3::el5"

#args <- simulateArgs(cmdLine)
args <- simulateArgs(cmdLine22)

# Print the args to verify
print(args)

## Import all.txt files in the current directory
#args = commandArgs(TRUE)
infile = args[1]  # data.csv
outDir = args[3]
clusterfile = args[4]  # cluster file
donorC = args[5]  # initial cell type
species = args[6]
orgname = args[7]
setwd("~/clean/transsynw")
## Load necessary libraries
source('code/SynergisticCore_library_noExpCutoff_pioneers_synergyBG.R') # 5 + all 3-4-5 nonJSDs
source('code/scJSD.R')
library(Rcpp)
sourceCpp('code/SynergisticCore.cpp')
sourceCpp('code/JSD.cpp')
require(gtools)
library(Matrix)
library(tibble)
library(dplyr)
library(stringr)
library(purrr)
library(data.table)
target.classes = str_split(args[2], "::", simplify = TRUE) %>% as.vector
merge=target.classes[length(target.classes)]
target.classes = target.classes[-length(target.classes)][-1]
orgname = str_split(args[7], "::", simplify = TRUE) %>% as.vector

print(target.classes)
print(merge)
print (species)
print(orgname)

dir.exists(outDir)  # Should return TRUE if directory exists
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)



## Load start pop
if(! file.exists((paste0("startpop/",donorC)))){
  geneexp=orgname[1]
  clust=orgname[2]
  donor_org=orgname[3]
  inipop = read.delim(donorC,sep="\t",header=TRUE,stringsAsFactors=FALSE,row.names=1)
  inipop = log2(inipop+1)
  donorE = rowMeans(inipop)
  #save(donorE, file='initial_gene_exprefssion.Robj')
} else {
  donor_org=donorC
  geneexp=orgname[1]
  clust=orgname[2]
  load((paste0("startpop/",donorC)))
}

## Load data
ftype = system(paste0("file -b --mime-type ",infile," | sed 's|/.*||'"),intern=TRUE)
if(ftype=='text'){
  M_all <- fread(file=args[1])
  M_all <- as.data.frame(M_all)
  rownames(M_all) <- M_all[[1]]
  M_all <- M_all[, -1]
} else {
  M_all = readRDS(infile)
}
rownames(M_all) = toupper(rownames(M_all))

# TF subsetting
tfs = read.csv(paste0("code/TF_",species,".txt"),header=FALSE,sep="\t",quote="")
M = M_all[toupper(rownames(M_all)) %in% toupper(tfs[,1]),]
print(dim(M))

## Assign cluster IDs
#map = read.delim(clusterfile,sep=",",header=FALSE,stringsAsFactors=FALSE)
#above is wrong
map <- fread(file=args[4])
map <- as.data.frame(map)
#rownames(map) <- map[[1]]
#map <-map[nrow(map):1,]
#map <- map[, -1]
map[2, ] <-sapply(map[2, ], function(x) gsub(" ", "_", x))

#merge_neurons <- TRUE
merge_neurons <- FALSE
#merge_neurons_except_target <- FALSE
merge_neurons_except_target <- TRUE
if (merge_neurons == T) {
  replace_categories <- c("L6_IT", "L6_CT", "L5_IT", "L2/3_IT", "L5/6_NP", "L5_ET", 
                          "Sst", "Vip", "Pvalb", "Lamp5", "L6b", "Sncg")
  
  # Sample vector `map[2,]` where you want to apply the changes
  map[2, ] <- sapply(map[2, ], function(x) {
    if (x %in% replace_categories) {
      x <- "L5_ET"
    } else {
      x <- x
    }
    x
  })
}
if (merge_neurons_except_target == T) {
  replace_categories <- c("L6_IT", "L6_CT", "L5_IT", "L2/3_IT", "L5/6_NP", 
                          "Sst", "Vip", "Pvalb", "Lamp5", "L6b", "Sncg")
  
  # Sample vector `map[2,]` where you want to apply the changes
  map[2, ] <- sapply(map[2, ], function(x) {
    if (x %in% replace_categories) {
      x <- "otherN"
    } else {
      x <- x
    }
    x
  })
}


m = match(colnames(M), as.character(map[1,]))
if(sum(is.na(m))==0){
  colnames(M) = map[2,m]
}else{
  stop("cluster-unassigned cells exist")
}
# merge subpopulations
if(merge=="yes"){
  colnames(M)[colnames(M) %in% target.classes] = paste0(target.classes,collapse="_")
  target.classes = paste0(target.classes,collapse="_")
}
unwanted_categories <- c("UNKNOWN", "DOUBLET", "LOW_QUALITY")
M <- M[, !toupper(colnames(M)) %in% unwanted_categories]
## All classes
classes = unique(colnames(M)); print(classes)

neurons <- M[, toupper(colnames(M)) %in% 'L5_ET']

M <- M[, !toupper(colnames(M)) %in% unwanted_categories]

###COMMENT AWAY IF NOT VS NEURON ONLY
# keep_categories <- c("L6_IT", "L6_CT", "L5_IT", "L2/3_IT", "L5/6_NP", "L5_ET",
#                     "Sst", "Vip", "Pvalb", "Lamp5", "L6b", "Sncg", "Astro")
# 
# # Subset the matrix to keep only columns whose names are in keep_categories
# M <- M[, colnames(M) %in% keep_categories]
###COMMENT AWAY IF NOT VS NEURON ONLY

colnames(M) <- sub("\\.\\d+$", "", colnames(M))


print(head(colnames(M)))



## Identify most synergistic cores.
for(cla in target.classes){
  if(cla != 'undefined'){ 
    print(cla)
    if (merge_neurons ==T){
      print("L5_ET means neurons in general")}
    computeMostSynergisticCores(data=M, classes=cla, infile=infile, donorC=donorC, outDir=outDir)
  }
}

## output
fil2 = paste0(outDir,"/analysis_summary.txt")
write.table(paste("Starting cell type:",donor_org, sep=" "), file=fil2,sep="\t",row.names=FALSE,quote=FALSE, col.names=FALSE,append=TRUE)
write.table(paste("Gene expression matrix:",geneexp, sep=" "), file=fil2,sep="\t",row.names=FALSE,quote=FALSE, col.names=FALSE,append=TRUE)
write.table(paste("Cluster file:",clust, sep=" "), file=fil2,sep="\t",row.names=FALSE,quote=FALSE, col.names=FALSE,append=TRUE)
write.table(paste("Species:",species, sep=" "), file=fil2,sep="\t",row.names=FALSE,quote=FALSE, col.names=FALSE,append=TRUE)
write.table(paste("Analysed subpopulations:",paste0(target.classes,collapse = ", "), sep=" "), file=fil2,sep="\t",row.names=FALSE,quote=FALSE, col.names=FALSE,append=TRUE)

core_files = list.files(path=outDir, pattern = "core",full.names = TRUE)
cores = map_dfr(core_files, function(x) {
  read.delim(x, sep = ',', stringsAsFactors = F)
})
write.table(cores, paste0(outDir,"/cores.tsv"), sep="\t", row.names = F)
unlink(core_files)

## Markers
gens = read.csv('code/all_marker_genes.txt',sep="\t",quote="",header=FALSE)
if(sum(toupper(rownames(M_all)) %in% toupper(gens[,1])) < 1000 ){stop("gene expression matrix does not contain enough marker genes")}

print(">Marker computation")
M = M_all[toupper(rownames(M_all)) %in% toupper(gens[,1]), ]
print(dim(M))
M = M[rowSums(M)!=0,,drop=FALSE]

# sizeNormalization
M = sweep(M,2,colSums(M),`/`)
M[is.na(M)] = 0
M = M[rowSums(M)!=0,,drop=FALSE]


indices <- match(colnames(M), map[1,])
colnames(M) <- map[2, indices]


JSD_rank_threshold = 10
for(cla in target.classes){
  if(cla != 'undefined'){
    print("Running C++ JSD...")
    JSD_value = scJSD(M, cla)
    sJSD_expr_score = scJSD_score(JSD_value)
    temp = sJSD_expr_score[1:JSD_rank_threshold]
    temp = cbind('Gene'=names(temp),'scJSD'=round(temp,2))
    fil = paste0(outDir,"/markers_",cla,".txt")
    write.table(temp, file=fil,sep=",",row.names=FALSE,quote=FALSE, col.names=TRUE)
  }	
}

# output
marker_files = list.files(path=outDir, pattern = "markers", full.names = TRUE)
markers = map_dfr(marker_files, function(x) {
  temp_df = read.delim(x, sep = ',', stringsAsFactors = F)
  pop = str_match(x, "markers_(.+)\\.txt")[2]
  temp_df = temp_df %>% mutate(Subpopulation = pop)
})
write.table(markers, paste0(outDir,"/markers.tsv"), sep="\t", row.names = F)
unlink(marker_files)


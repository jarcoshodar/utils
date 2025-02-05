#!/bin/bash                                       

#because chip seq atlases are too big to keep all in memory outside workstations/hpc, it makes sense to subset to what your query is known to be
#in our case, the known tfs in a differentially gene expressed profile (as in wrapper folders)
#needs modification if your dge has different column order/names, or you prefer other sginificance threshold
function extract_de_tfs() {
    #usage: extract_de_tfs <tf_list_file> <dge_results_file> <output_file>
    
    
    if [ "$#" -ne 3 ]; then
        echo "Error: Incorrect number of arguments"
        echo "Usage: extract_de_tfs <tf_list_file> <dge_results_file> <output_file>"
        echo "Example: extract_de_tfs Mus_musculus_TF_LIST.txt dge_results.csv output_tfs.txt"
        return 1
    fi

    local tf_file="$1"
    local dge_file="$2"
    local output_file="$3"

    if [ ! -f "$tf_file" ] || [ ! -f "$dge_file" ]; then
        echo "Error: Input file(s) not found"
        return 1
    fi

    if [ ! -r "$tf_file" ] || [ ! -r "$dge_file" ]; then
        echo "Error: Input file(s) not readable"
        return 1
    fi

    #detect separator for DGE file
    local separator=$(grep -q "," "$dge_file" && echo "," || echo "\t")

    #get significant DE TFs
    #modify this if your dge has a different format
    #currently indexes on column number and 0.05 on 7 = p val adj and gene symbols on 2nd
    awk -v sep="$separator" -F"[,\t]" '
        NR==FNR && FNR>1 {tf[$2]=1; next}
        FNR>1 && $7<0.05 && $2 in tf {print $2}
    ' "$tf_file" "$dge_file" > "$output_file"
    if [ $? -eq 0 ] && [ -s "$output_file" ]; then
        echo "Successfully created $output_file with $(wc -l < "$output_file") transcription factors"
        return 0
    else
        echo "Error: Problem creating output file"
        return 1
    fi
}

#function to create a subset bed file based on a gene list
#can be be does not have to be as above
#that being said, doubtful it works if your .bed isn't chip seq atlas
function subset_bed_by_genes() {
    #usage: subset_bed_by_genes <gene_list> <input_bed> <output_bed>
    
    if [ "$#" -ne 3 ]; then
        echo "Error: Incorrect number of arguments"
        echo "Usage: subset_bed_by_genes <gene_list> <input_bed> <output_bed>"
        echo "Example: subset_bed_by_genes DE_TFs.txt full.bed subset.bed"
        return 1
    fi

    local gene_list="$1"
    local input_bed="$2"
    local output_bed="$3"

    if [ ! -f "$gene_list" ] || [ ! -f "$input_bed" ]; then
        echo "Error: Input file(s) not found"
        return 1
    fi

    if [ ! -r "$gene_list" ] || [ ! -r "$input_bed" ]; then
        echo "Error: Input file(s) not readable"
        return 1
    fi

    #create pattern and subset bed file
    local pattern=$(awk '{print "Name="$1"%20\\("}' "$gene_list" | paste -sd'|' -)
    grep -E "$pattern" "$input_bed" > "$output_bed"
    if [ $? -eq 0 ] && [ -s "$output_bed" ]; then
        echo "Successfully created $output_bed with $(wc -l < "$output_bed") regions"
        return 0
    else
        echo "Error: Problem creating output bed file"
        return 1
    fi
}

# Usage examples:
# First source this file
# or add the functions to your bashrc
# source tf_functions.sh
#
# Then you can use either function independently:
# extract_de_tfs Mus_musculus_TF_LIST.txt dge_results.csv DE_TFs.txt
# subset_bed_by_genes DE_TFs.txt full.bed subset.bed
#
# Or pipeline them together:
# extract_de_tfs Mus_musculus_TF_LIST.txt dge_results.csv DE_TFs.txt && \
# subset_bed_by_genes DE_TFs.txt Oth.ALL.50.AllAg.AllCell.bed subset_TF_regions.bed

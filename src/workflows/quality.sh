#!/bin/bash

# ==============================================================================
#? Sample quality control of PANC-1 RNA-Seq expression data
# ==============================================================================
 
# --- General settings and variables -------------------------------------------
source ./src/bash_commons.sh

in_path="./data/in/"
out_path="./data/out/expression/"
expression_matrix="pH_CountMatrix_genes_TPM.tsv"
out_sub="QC"

# --- The pipeline starts here -------------------------------------------------
echo -e "\n${mag}STARTING DATA QUALITY CONTROL${end}"

# Check x.FASTQ is properly installed
if ! which x.fastq > /dev/null 2>&1; then
    echo "Cannot find x.FASTQ..."
    echo -e "${red}\nPIPELINE FAILED${end}\n"
    exit 1
else
    echo "x.FASTQ found!"
fi

# Check if target directory exists
if [ ! -d "$out_path" ]; then
  mkdir -p "$out_path"
fi

# This cp is because qcfastq, by default, saves output in the same folder
# used to locate the input files 
cp "${in_path}/${expression_matrix}" "${out_path}/${expression_matrix}"
qcfastq --workflow --tool=PCA --suffix="TPM.tsv" --out="$out_sub" "$out_path"

# Do the cleaning
rm "${out_path}/${expression_matrix}"
log_file="$(ls "${out_path}"/Z_QC_PCA*.log)"
mv "${log_file}" "${out_path}/${out_sub}/$(basename ${log_file})"

# --- The pipeline ends here ---------------------------------------------------
if [[ $? -eq 0 ]]; then
    echo -e "${mag}PIPELINE COMPLETED SUCCESSFULLY${end}"
else
    echo -e "\n${red}PIPELINE FAILED${end}"
fi

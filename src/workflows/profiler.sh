#!/bin/bash

# ==============================================================================
#? Channelome absolute expression analysis in PANC-1 cell line
# ==============================================================================
 
# --- General settings and variables -------------------------------------------
source ./src/bash_commons.sh

in_path="./data/in/"
out_path="./data/out/expression/"
GOIs="./data/in/pore_set.csv"
expression_matrix="pH_CountMatrix_genes_TPM.tsv"

# --- The pipeline starts here -------------------------------------------------
echo -e "\n${mag}STARTING PANC-1 PROFILING${end}"

# Check if target directory exists
if [ ! -d "$out_path" ]; then
  mkdir -p "$out_path"
fi

# Channelome absolute expression profiling in PANC-1 cell line
echo "Running panc_profiler.R ..."
Rscript --vanilla "./src/panc_profiler.R" \
    "${in_path}/${expression_matrix}" \
    "$GOIs" \
    "$out_path"

# --- The pipeline ends here ---------------------------------------------------
if [[ $? -eq 0 ]]; then
    echo -e "${mag}PIPELINE COMPLETED SUCCESSFULLY${end}"
else
    echo -e "\n${red}PIPELINE FAILED${end}"
fi

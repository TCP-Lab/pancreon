#!/bin/bash

# ==============================================================================
#? Whole-cell current I-V analysis in PANC-1 cell line
# ==============================================================================
 
# --- General settings and variables -------------------------------------------
source ./src/bash_commons.sh

in_path="./data/in/patch"
out_path="./data/out/"

# --- The pipeline starts here -------------------------------------------------
echo -e "\n${mag}STARTING WHOLE-CELL CURRENT ANALYSIS${end}"
echo -e "${mag} ====================================${end}"

# Whole-cell current IV analysis in PANC-1 cell line
echo "Running compare_ramps.R ..."
Rscript --vanilla "./src/compare_ramps.R" \
    "${in_path}" \
    "$out_path"

# --- The pipeline ends here ---------------------------------------------------
if [[ $? -eq 0 ]]; then
    echo -e "\n${mag}===============================${end}"
    echo -e "${mag}PIPELINE COMPLETED SUCCESSFULLY\n${end}"
else
    echo -e "${red}\nPIPELINE FAILED${end}\n"
fi

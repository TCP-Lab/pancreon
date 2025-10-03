#!/bin/bash

# ==============================================================================
#? Channelome geneset build by querying the MTP-DB (v1.25.24)
#?
#? Transportome sub-setting is done according to the following criteria:
#? - ICs: all (436)
#? - AQPs: all (14)
#? - SLCs: none
#? - Pumps: none
#? - ABCs: none
# ==============================================================================

# --- General settings and variables -------------------------------------------
source ./src/bash_commons.sh

db_path="./data/MTPDB.sqlite"
out_dir="./data/in"

# --- The pipeline starts here -------------------------------------------------
echo -e "\n${mag}STARTING GOI LIST ASSEMBLY${end}"

# --- Extract the archive ------------------------------------------------------
_extract_mtpdb "$db_path"

# --- Make the gene-set --------------------------------------------------------
printf "Making the geneset...\n"
Rscript --vanilla "./src/make_geneset.R" \
    "$db_path" \
    "$out_dir"

# --- The pipeline ends here ---------------------------------------------------
if [[ $? -eq 0 ]]; then
    echo -e "\n${mag}PIPELINE COMPLETED SUCCESSFULLY${end}"
else
    echo -e "\n${red}PIPELINE FAILED${end}"
fi

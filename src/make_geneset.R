#!/usr/bin/env Rscript

# --- Packages -----------------------------------------------------------------

#library(dplyr, warn.conflicts = FALSE)
library(DBI)
library(RSQLite)

# --- Functions ----------------------------------------------------------------

# --- Input Parsing ------------------------------------------------------------

# Extract command-line arguments
db_path <- commandArgs(trailingOnly = TRUE)[1]
out_dir <- commandArgs(trailingOnly = TRUE)[2]

# # Live debug (from the project root directory)
# db_path <- "./data/MTPDB.sqlite"
# out_dir <- "./data/in"

# if (Sys.info()["sysname"] == "Windows") {
#   # Can't access the DB directly on WSL when running from Windows...
#   db_path <- file.path(Sys.getenv("USERPROFILE"), "Desktop", "MTPDB.sqlite")
#   file.copy(from = "./data/MTPDB.sqlite", to = db_path)
# }

# Connect to the MTP-DB
connection <- dbConnect(SQLite(), dbname = db_path)

# Pores (Ion Channels + AQPs)
query_pores <-
  "SELECT DISTINCT
    	channels.ensg,
    	gene_names.hugo_gene_symbol,
    	gene_names.hugo_gene_name
    FROM
    	channels JOIN gene_names ON channels.ensg = gene_names.ensg
    
    UNION
    
    SELECT DISTINCT
    	aquaporins.ensg,
    	gene_names.hugo_gene_symbol,
    	gene_names.hugo_gene_name
    FROM
    	aquaporins JOIN gene_names ON aquaporins.ensg = gene_names.ensg
    
    ORDER BY gene_names.hugo_gene_symbol"

# Make the calls
pores <- dbGetQuery(connection, query_pores)

# Extract Gene Symbols, combine, and save
pores$hugo_gene_symbol |>
  unique() |> na.omit() |>
  write.table(sep = ",",
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE,
              file = file.path(out_dir, "pore_set.csv"))

# Disconnect from the MTP-DB
dbDisconnect(connection)






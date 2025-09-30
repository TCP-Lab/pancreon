# --- Packages -----------------------------------------------------------------

library(ggplot2)
library(r4tcpl)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
#library(DBI)
#library(RSQLite)
# library(httr)



# Function loading
source("./src/util_functions.R")

gois_file <- "./data/out/pore_set.csv"
expression_file <- "./data/in/pH_CountMatrix_genes_TPM.tsv"
out_dir <- "./data/out"

# Load TPMs
expression_file |> read.delim(header = TRUE) -> expression_matrix

# Load the list of GOIs
gois_file |> read.delim(header = FALSE) -> gois


# Subsetting
#expression_matrix[expression_matrix$SYMBOL %in% unlist(gois),] |> dim()
expression_matrix |> filter(SYMBOL %in% unlist(gois)) -> gois_expression


# Take the log2
gois_expression[,c(5:16)] <- log2(gois_expression[,c(5:16)] + 1)


# Control condition

gois_expression |> select(1:3, ends_with("Ctrl")) |> mutate(
  Mean = rowMeans(across(ends_with("Ctrl")), na.rm = TRUE),
  Std_Dev = apply(across(ends_with("Ctrl")), 1, sd, na.rm = TRUE)) |> filter(Mean >= 1) -> ctrl_mean


plot_barChart(ctrl_mean, y_limit = 10, border = FALSE, thr = 1, out_folder = out_dir)


  
  
# Selected condition

gois_expression |> select(1:3, ends_with("Select")) |> mutate(
  Mean = rowMeans(across(ends_with("Select")), na.rm = TRUE),
  Std_Dev = apply(across(ends_with("Select")), 1, sd, na.rm = TRUE)) |> filter(Mean >= 1) -> select_mean


plot_barChart(select_mean, y_limit = 10, border = FALSE, thr = 1, out_folder = out_dir)
  


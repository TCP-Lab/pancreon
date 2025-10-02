#!/usr/bin/env Rscript

# --- Packages -----------------------------------------------------------------

library(ggplot2)
#library(r4tcpl)
library(dplyr, warn.conflicts = FALSE)

# Function loading
source("./src/util_functions.R")

# --- Input Parsing ------------------------------------------------------------

# Extract command-line arguments.
expression_file <- commandArgs(trailingOnly = TRUE)[1]
gois_file <- commandArgs(trailingOnly = TRUE)[2]
out_dir <- commandArgs(trailingOnly = TRUE)[3]

# # Interactive debug (from the project root directory)
# expression_file <- "./data/in/pH_CountMatrix_genes_TPM.tsv"
# gois_file <- "./data/in/pore_set.csv"
# out_dir <- "./data/out"

# --- Data Loading and Pre-processing ------------------------------------------

# Load Expression Data (in TPMs)
expression_file |> read.delim(header = TRUE) -> expression_matrix

# Load the list of GOIs
gois_file |> read.delim(header = FALSE) -> gois

# Subset Expression Matrix
expression_matrix |> filter(SYMBOL %in% unlist(gois)) -> gois_expression

# Take the log2(TPM+1), overwriting the original columns
group_labels <- c("Ctrl", "Acute", "Select")
regex <- paste0("\\.(", paste(group_labels, collapse = "|"), ")$")
gois_expression |> mutate(across(matches(regex), \(x){log2(x+1)})) -> gois_expression

# --- Make Statistics ----------------------------------------------------------

# Parameter Settings
y_limit <- 10
border <- FALSE
thr <- 1

for (condition in group_labels) {
    # Compute basic descriptive stats (Mean and SD)
    gois_expression |> select(1:3, ends_with(condition)) |>
        mutate(
            Mean = rowMeans(across(ends_with(condition)), na.rm = TRUE),
            Std_Dev = apply(across(ends_with(condition)), 1, sd, na.rm = TRUE)
        ) |> filter(Mean >= 1) -> goi_stats
    # Save Stat Table
    write.csv(goi_stats,
              file.path(out_dir, paste0("channelome_stats_", condition, ".csv")),
              row.names = FALSE)
    
    # Print Plots
    plot_barChart(goi_stats, condition, y_limit, border, thr) -> barChart
    # Save the Chart
    r4tcpl::savePlots(
        \(){print(barChart)},
        width_px = 2000,
        figure_Name = paste0("channelome_chart_", condition),
        figure_Folder = out_dir)
}




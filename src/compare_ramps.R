#!/usr/bin/env Rscript

# --- Packages -----------------------------------------------------------------

library(readABF)
#library(r4tcpl)

# Function loading
source("./src/util_functions.R")

# --- Input Parsing ------------------------------------------------------------

# Extract command-line arguments.
in_dir <- commandArgs(trailingOnly = TRUE)[1]
out_dir <- commandArgs(trailingOnly = TRUE)[2]

# # Interactive debug (from the project root directory)
# in_dir <- "./data/in/patch"
# out_dir <- "./data/out"

# Create a named vector with the ending sample of each protocol phase
# (i.e., duration of the Holding Potential in time samples, ad so on...)
proto <- c(holding = 32, pre = 533, ramp = 1532)

# --- Data Loading -------------------------------------------------------------

# List all subfolders (each one = experiment)
subfolders <- list.dirs(in_dir, full.names = TRUE, recursive = FALSE)

# Loop over experiments (subfolders)
for (exp_path in subfolders) {
  
  # List all ABF files in this experiment folder
  abf_files <- list.files(exp_path, pattern = "\\.abf$", full.names = TRUE)
  if (length(abf_files) == 0) {
    warning(paste("No ABF files found in", exp_path, "input directory."))
    next
  }
  
  # Get the experiment ID
  exp_id <- basename(exp_path)
  
  # Initialize an empty list (of data frames)
  ramps <- list()
  # Loop over ABF files and store them in the 'ramps' list
  for (f in abf_files) {
    
    basename(f) |>
      sub(paste0(exp_id, "_"), "", x=_, fixed = TRUE) |>
      sub(".abf", "", x=_, ignore.case = FALSE, fixed = TRUE) -> label
    
    ramps[[label]] <- as.data.frame(readABF(f))
  }
  
  # Check the compatibility of the ramps before building the reduced dataframe
  ramps |> compatibility_test(t_heading = "Time [s]", t_err = 1e-5,
                              v_heading = "10_Vm [mV]", v_err = 1)
  
  # Collapse list into a single data frame with one time and voltage columns
  df <- data.frame(Time = ramps[[1]]$`Time [s]`, Vm = ramps[[1]]$`10_Vm [mV]`)
  for (label in names(ramps)) {
    df[[label]] <- ramps[[label]]$`Im_scaled [pA]`
  }
  
  # --- Pre-processing ---------------------------------------------------------
  
  # Plot the Voltage Protocol
  r4tcpl::savePlots(\(){
    plot_voltage_protocol(df, exp_id)},
    width_px = 1000,
    figure_Name = "Voltage_Protocol",
    figure_Folder = file.path(out_dir, exp_id))
  
  # Plot Full Ramps in Time
  r4tcpl::savePlots(\(){
    plot_full_ramps(df, proto, exp_id)},
    width_px = 1000,
    figure_Name = "All_Ramps",
    figure_Folder = file.path(out_dir, exp_id))
  
  # Plot Leakage under holding potential
  r4tcpl::savePlots(\(){
    plot_holding_leakage(df, proto, exp_id)},
    width_px = 1000,
    figure_Name = "Holding_Leakage",
    figure_Folder = file.path(out_dir, exp_id))
  
  # --- I-V Analysis -----------------------------------------------------------
  
  # Load the contrasts of interest
  list.files(exp_path,
             pattern = "^comp.*\\.txt$",
             full.names = TRUE,
             ignore.case = TRUE) -> comp_files
  
  if (length(comp_files) == 1) {
    comp_files |> read.delim(header = FALSE) |> unlist() -> comparisons
  } else if (length(comp_files) == 0) {
    cat("WARNING: comparison definition unavailable... skip this eperiment")
    next
  } else if (length(comp_files) > 1) {
    cat("WARNING: multiple comparison definitions available: taking one...")
    comp_files[1] |> read.delim(header = FALSE) |> unlist() -> comparisons
  }
  
  # Make all the planned comparisons
  for (comp in comparisons) {
    # Parse the contrast
    comp |> strsplit("-") |> unlist() |> {\(x)x[1]}() -> cond
    comp |> strsplit("-") |> unlist() |> {\(x)x[2]}() -> ctrl
    
    # Plot the I-V curves for the conditions of interest
    r4tcpl::savePlots(\(){
      plot_diff_IV(df, proto, cond, ctrl, exp_id)},
      width_px = 1000,
      figure_Name = paste0("IV_", cond, "-", ctrl),
      figure_Folder = file.path(out_dir, exp_id))
  }
}


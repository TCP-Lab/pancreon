

# My colors
nero <- "#05072b"
blu <- "#1146cf"
rosso <- "#de1620"

# 'gois_stats' is supposed to be a data.frame with:
#   - descriptive stats of gene expression in 'Mean' and 'Std_Dev' columns;
#   - a gene annotation columns featuring (at least) the 'SYMBOL' key;
plot_barChart <- function(gois_stats,
                          data_label,
                          y_limit,
                          border = FALSE,
                          thr)
{
  # Colors
  line_color <- "gray17"
  err_color <- "gray17"
  # Prepare the Frame
  gg_frame <-
    ggplot(data = gois_stats,
           aes(x = SYMBOL, y = Mean, fill = SYMBOL)) +
    theme_bw(base_size = 15, base_rect_size = 1.5) +
    theme(axis.text.x = element_text(size = 10, angle = 90,
                                     vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.position = "none") +
    scale_y_continuous(limits=c(0, y_limit),
                       expand = c(0.01, 0.02),
                       breaks = seq(0, y_limit, 0.5)) +
    xlab("Genes of Interest") +
    ylab(substitute(log[2]*(x+1), list(x = "TPM"))) +
    ggtitle(label = paste("Expressed Channelome in", data_label, "condition"))
  
  # Draw the Bars
  if (border) {
    gg_bars <- geom_bar(stat = "identity", width = 0.7,
                        color = line_color, linewidth = 0.1)
  } else {
    gg_bars <- geom_bar(stat = "identity", width = 0.75)
  }
  gg_errorbar <- gg_frame + gg_bars +
    geom_errorbar(aes(ymin = Mean - Std_Dev,
                      ymax = Mean + Std_Dev),
                  linewidth = 1.0, width = 0.5, color = err_color)
  # Add the Expression Threshold
  gg_thr <- gg_errorbar +
    geom_hline(yintercept = thr,
               linetype = "dashed",
               color = line_color,
               linewidth = 1)
  
  return(gg_thr)
}

# This function returns expanded 'ylim' for better plotting, based on the
# overall range of values of the input data frame and the 'increase' parameter:
# 0.1 == +10% increase
expand <- function(df, increase)
{
  # Total y span
  df |> unlist() |> range() -> y_lims
  delta <- y_lims[2] - y_lims[1]
  
  extra <- delta*increase
  
  return(c(y_lims[1]-extra, y_lims[2]+extra))
}

# Plot the Voltage Protocol (for a dataframe with 'Time' and 'Vm' columns)
plot_voltage_protocol <- function(df, exp_id)
{
  # Save current settings and fix margins
  old_par <- par(no.readonly = TRUE)
  par(mar = c(5.1, 5.5, 4.1, 2.5)) # c(bottom, left, top, right)
  # Make the plot
  plot(df$Time * 1e3, df$Vm,
       type = "l",
       lty = 1,
       lwd = 2,
       main = "Voltage Protocol",
       xlab = "time (ms)",
       ylab = substitute(V[m]*" "*(mV)),
       cex.main = 2,      # Main title size
       cex.lab = 1.5,     # Axis labels size
       cex.axis = 1.2)    # Axis tick label size
  # Add Subtitle
  mtext(paste("Experiment:", exp_id),
        side = 3,
        line = -1.5,
        cex = 1.2)
  # Restore original settings
  par(old_par) |> suppressWarnings()
}

# Plot all the current ramps in a dataframe, as a function of the 'Time' column.
# 'proto' is a named vector with the ending sample of each protocol phase, in
# particular a 'pre' and 'ramp' value, representing the beginning and the end of
# the voltage sweep, respectively
plot_full_ramps <- function(df, proto, exp_id)
{
  # Only currents
  currents <- df[,-c(1,2)]
  
  # Magnify the ramp phase, with a 10% margin (cut the capacitive transients)
  expand(currents[proto["pre"]:proto["ramp"],], increase = 0.1) -> y_lims
  
  # Save current settings and fix margins
  old_par <- par(no.readonly = TRUE)
  par(mar = c(5.1, 5.5, 4.1, 2.5)) # c(bottom, left, top, right)
  # Make the plot
  matplot(df$Time * 1e3, currents,
          type = "l",
          lty = 1,
          lwd = 2,
          col = rainbow(ncol(currents)),
          ylim = y_lims,
          main = "Current Ramps",
          xlab = "time (ms)",
          ylab = substitute(I[m]*" "*(pA)),
          cex.main = 2,      # Main title size
          cex.lab = 1.5,     # Axis labels size
          cex.axis = 1.2)    # Axis tick label size
  # Add Subtitle
  mtext(paste("Experiment:", exp_id),
        side = 3,
        line = -1.5,
        cex = 1.2)
  legend("topright",
         legend = colnames(currents),
         y.intersp = 1.5,
         col = rainbow(ncol(currents)),
         lty = 1,
         lwd = 2,
         pch = 19,
         bty = "n")
  # Restore original settings
  par(old_par) |> suppressWarnings()
}

# Plot the beginning of each ramp to analyze the average leakage current as a
# function of the 'Time' column.
# 'proto' is a named vector with the ending sample of each protocol phase, in
# particular a 'holding' value, representing the end of the holding voltage.
plot_holding_leakage <- function(df, proto, exp_id)
{
  # Only leakage currents
  leak <- df[1:proto["holding"],-c(1,2)]
  
  # Magnify the leakage, with a 10% margin
  expand(leak, increase = 0.1) -> y_lims
  
  # Make the plot
  # Save current settings and fix margins (to make room for the legend)
  old_par <- par(no.readonly = TRUE)
  par(mar = c(5.1, 5.5, 4.1, 12)) # c(bottom, left, top, right)
  matplot(df$Time[1:proto["holding"]] * 1e3, leak,
          type = "l",
          lty = 1,
          lwd = 2,
          col = rainbow(ncol(leak)),
          ylim = y_lims,
          main = "Leakage @ Holding Potential",
          xlab = "time (ms)",
          ylab = substitute(I[m]*" "*(pA)),
          cex.main = 2,      # Main title size
          cex.lab = 1.5,     # Axis labels size
          cex.axis = 1.2)    # Axis tick label size
  # Add Subtitle
  mtext(paste("Experiment:", exp_id),
        side = 3,
        line = -1.5,
        cex = 1.2)
  legend("right",
         legend = paste0(colnames(leak), ": ", round(colMeans(leak),1), " pA"),
         text.font = 1,
         y.intersp = 2,
         col = rainbow(ncol(leak)),
         lty = 1,
         lwd = 2,
         pch = 19,
         inset = -0.2,
         xpd = TRUE,
         bty = "n")
  # Restore original settings
  par(old_par) |> suppressWarnings()
}

# Plot the I-V curves of two conditions of interest, together with their
# difference (a 'Vm' column is expected in the input dataframe 'df').
# 'proto' is a named vector with the ending sample of each protocol phase, in
# particular a 'pre' and 'ramp' value, representing the beginning and the end of
# the voltage sweep, respectively.
# 'cond' and 'ctrl' are the label of the conditions of interest used to select
# columns from the input dataframe 'df'.
plot_diff_IV <- function(df, proto, cond, ctrl, exp_id)
{
  # Cut across the ramp 
  df_cut <- df[proto["pre"]:proto["ramp"],]
  
  # Estimate reversal potential
  rev_ctrl <- reversal_potential(df_cut$Vm, df_cut[[ctrl]], 5, exp_id)
  rev_cond <- reversal_potential(df_cut$Vm, df_cut[[cond]], 5, exp_id)
  
  # Compute the difference I-V
  difference <- df_cut[[cond]] - df_cut[[ctrl]]
  
  # Save current settings and fix margins
  old_par <- par(no.readonly = TRUE)
  par(mar = c(5.1, 5.5, 4.1, 2.5)) # c(bottom, left, top, right)
  # Make the plot
  matplot(df_cut$Vm,
          cbind(df_cut[[ctrl]], df_cut[[cond]], difference),
          type = "l",
          lty = 1,
          lwd = c(2, 2, 3),
          col = c(nero, blu, rosso),
          main = paste("I-V Difference [", cond, "\u2212", ctrl, "]"), # \u2212 for minus sign
          xlab = substitute(V[m]*" "*(mV)),
          ylab = substitute(I[m]*" "*(pA)),
          cex.main = 2,      # Main title size
          cex.lab = 1.5,     # Axis labels size
          cex.axis = 1.2)    # Axis tick label size
  # Add Reversal Potentials
  usr <- par("usr") # Get plot limits, where usr[4] is the upper y-limit
  text(rev_ctrl, usr[4]/5,
       substitute(V[rev] == x*" "*mV, list(x = round(rev_ctrl,1))),
       col = nero, cex = 1)
  segments(rev_ctrl, usr[4]/50,
           rev_ctrl, usr[4]/5 - usr[4]/35,
           col = nero, lwd = 1, lty = 3)
  text(rev_cond, usr[4]/3,
       substitute(V[rev] == x*" "*mV, list(x = round(rev_cond,1))),
       col = blu, cex = 1)
  segments(rev_cond, usr[4]/50,
           rev_cond, usr[4]/3 - usr[4]/30,
           col = blu, lwd = 1, lty = 3)
  # Add Subtitle
  mtext(paste("Experiment:", exp_id),
        side = 3,
        line = -1.5,
        cex = 1.2)
  lines(df_cut$Vm,
        rep(0, length(df_cut$Time)),
        col = "black",
        lwd = 1,
        lty = 5)
  legend("topleft",
         inset = 0.03,
         legend = c(ctrl, cond, "delta"),
         y.intersp = 1.5,
         col = c(nero, blu, rosso),
         lty = 1,
         lwd = 2,
         pch = 19,
         bty = "n")
  # Restore original settings
  par(old_par) |> suppressWarnings()
}

# Check if two vectors are equal to within an user-defined arbitrary error
equal_to_within_err <- function(x, y, err = 1e-1)
{
  max(abs(x-y)) < err
}

# Test if all the ramps within a list have the same voltage protocol
compatibility_test <- function(ramps,
                               t_heading = "Time [s]",
                               t_err = 1e-5,
                               v_heading = "10_Vm [mV]",
                               v_err = 1)
{
  err_msg <- "ERROR: Different voltage protocols detected!"
  
  # Check vector lengths
  ramps |> sapply(\(ramp){
    nrow(ramps[[1]]) == nrow(ramp)
  }) |> all() -> test_out
  if (!test_out) {cat(err_msg)}
  
  # Check time vectors
  ramps |> sapply(\(ramp){
    equal_to_within_err(ramps[[1]][[t_heading]],
                        ramp[[t_heading]],
                        err = t_err)
  }) |> all() -> test_out
  if (!test_out) {cat(err_msg)}
  
  # Check voltage vector
  ramps |> sapply(\(ramp){
    equal_to_within_err(ramps[[1]][[v_heading]],
                        ramp[[v_heading]],
                        err = v_err)
  }) |> all() -> test_out
  if (!test_out) {cat(err_msg)}
}

# Fit the I-V curve with an n-degree polynomial and find the real root to
# estimate the reversal potential
reversal_potential <- function(x, y, n = 5, exp_id = NULL)
{
  # Fit an n-degree polynomial
  model <- lm(y ~ poly(x, n, raw=TRUE))
  
  # For debug only
  if (FALSE) {
    # Predict at new points and plot
    x_new <- seq(x[1], x[length(x)], length.out = 100)
    y_pred <- predict(model, data.frame(x = x_new))
    lines(x_new, y_pred, col="red")
  }
  
  # Get coeeficients and find roots (zeroes)
  model |> coef() |> polyroot() -> roots
  # Keep only real roots
  real_roots <- Re(roots[abs(Im(roots)) < 1e-8])
  
  if (length(real_roots) > 1) {
    cat("WARNING: multiple roots detected in ", exp_id,
        "... only one is shown!\n", sep = "")
  }
  return(real_roots[1])
}


## --- Currently Unsued --------------------------------------------------------

# Regenerate annotation from scratch
add_annotation <- function(gene_matrix, OrgDb_key = "ENSEMBL") {
  # See columns(org.Hs.eg.db) or keytypes(org.Hs.eg.db) for a complete list of
  # all possible annotations.
  org_db <- org.Hs.eg.db::org.Hs.eg.db
  annots <- AnnotationDbi::select(org_db,
                                  keys = gene_matrix[,"IDs"],
                                  columns = c("SYMBOL", "GENENAME", "GENETYPE"),
                                  keytype = OrgDb_key)
  # Warning: 'select()' returned 1:many mapping between keys and columns
  # ========>
  # Collapse the duplicated entries in the ID column and concatenate the
  # (unique) values in the remaining columns using a comma as a separator to
  # prevent rows from being added in the following join step.
  #if (anyDuplicated(annots[,OrgDb_key])) {
  #  cat("\nWARNING:\n Multiple annotation entries corresponding to a single\n",
  #      OrgDb_key, "ID will be collapsed by a comma separator.\n")
  #  annots <- aggregate(. ~ get(OrgDb_key),
  #                      data = annots,
  #                      FUN = \(x)paste(unique(x), collapse = ","),
  #                      na.action = NULL)[,-1]
  #}
  colnames(annots)[colnames(annots) == OrgDb_key] <- "IDs" 
  gene_matrix <- merge(annots, gene_matrix, by = "IDs", all.y = TRUE)
}

# Correlation ScatterPlot Matrix (with fixed text size)
custom_pairs <- function(data_set, color = "gray15") {
  # Customize lower panel (correlation values)
  panel_cor <- function(x, y) {
    default_usr <- par("usr")
    on.exit(par(usr = default_usr))
    par(usr = c(0, 1, 0, 1))
    r <- round(cor(x, y), digits = 3)
    text(0.5, 0.5, r, cex = 5)
  }
  # Customize upper panel (scatter plots)
  panel_points <- function(x, y) {
    points(x, y, pch = 19, cex = 1, col = color)
  }
  # Create the plots
  par(cex.axis = 2)
  pairs(data_set,
        cex.labels = 4,
        font.labels = 4,
        lower.panel = panel_cor,
        upper.panel = panel_points)
}

profilePlot <- function(matrix_of_means, chart_type = "boxplot", thr = 1) {
    
    # Specify Data
    master_plot <-
        ggplot(matrix_of_means,
               aes(x = SYMBOL,
                   y = Mean,
                   group = SYMBOL))
    
    # Specify Graphic Elements
    if (chart_type == "boxplot") {
        master_plot <- master_plot +
            geom_boxplot(
                width = 0.7,
                fill = "mediumpurple1",
                color = "gray20",   # darker border
                alpha = 0.5) +
            geom_jitter(
                aes(color = Source),
                width = 0.2,
                size = 1,
                alpha = 0.6)
    } else if (chart_type == "95ci") {
        master_plot <- master_plot +
            stat_summary(
                fun.data  = mean_cl_normal,
                fun.args  = list(conf.int = 0.95),
                geom      = "errorbar",
                width     = 0.2,
                color     = "black",
                linewidth = 1) +
            # point at the mean
            stat_summary(
                fun   = mean,
                geom  = "point",
                size  = 3,
                color = "darkred")
    } else if (chart_type == "points") {
        master_plot <- master_plot +
            geom_point(
                aes(color = Source))
    } else if (chart_type == "lines") {
        master_plot <- master_plot +
            geom_line(
                aes(group = Source,
                    color = Source))
    } else {
        cat("\nUndefined 'chart_type'")
    }
    
    # Add a Threshold
    master_plot <- master_plot +
        geom_hline(yintercept = thr,
                   linetype = "dashed",
                   color = "gray17",
                   linewidth = 1)
    
    # Flip the plot
    master_plot <- master_plot +
        scale_x_discrete(limits = rev(levels(factor(matrix_of_means$SYMBOL)))) +
        coord_flip()
    
    # Specify Frame Attributes
    master_plot <- master_plot +
        theme_bw(base_size = 15, base_rect_size = 1.5) +
        theme(axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 10),
              axis.title = element_text(size = 18),
              legend.text = element_text(size = 16),
              legend.position = "inside",
              legend.position.inside = c(0.8, 0.1)) +
        xlab("Genes of Interest") +
        ylab(substitute(log[2]*(x+1), list(x = "TPM"))) +
        ggtitle(label = "Expression Profile Plot")
}

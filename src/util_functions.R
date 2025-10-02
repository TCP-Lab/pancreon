

# 'gois_stats' is supposed to be a data.frame with:
#   - descriptive stats on gene expression in 'Mean' and 'Std_Dev' columns;
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

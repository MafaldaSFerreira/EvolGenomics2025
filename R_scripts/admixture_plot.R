#!/usr/bin/env Rscript

# -------------------------
# LIBRARIES
# -------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

# -------------------------
# FUNCTION
# -------------------------
plot_admixture <- function(files, sample_order_file, metadata_file, output_png = "admixture_plots.png") {
  
  # Read inputs
  metadata <- read.table(metadata_file, header = TRUE)
  sample_order <- read.table(sample_order_file, header = FALSE)
  
  plots <- list()
  
  for (f in files) {
    # Read Q file
    Q_tbl <- read.table(f)
    
    # Add individuals from sample order
    Q_tbl$indv <- sample_order$V1
    
    # Merge with metadata
    Q_tbl_with_metadata <- merge(Q_tbl, metadata, by.x = "indv", by.y = "sample_id")
    
    # Reorder individuals by population
    Q_tbl_with_metadata <- Q_tbl_with_metadata[order(Q_tbl_with_metadata$population_code), ]
    Q_tbl_with_metadata$indv <- factor(Q_tbl_with_metadata$indv,
                                       levels = Q_tbl_with_metadata$indv)
    
    # Convert to long format
    Q_tbl_long <- Q_tbl_with_metadata %>%
      pivot_longer(cols = starts_with("V"),
                   names_to = "K",
                   values_to = "admix_proportion")
    
    # Make plot
    p <- ggplot(Q_tbl_long, aes(x = indv, y = admix_proportion, fill = K)) +
      geom_bar(position = "stack", stat = "identity") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            axis.title.x = element_blank()) +
      ggtitle(basename(f))
    
    plots[[f]] <- p
  }
  
  # Combine all plots in a grid
  combined_plot <- plot_grid(plotlist = plots, ncol = 2)  # adjust ncol as needed
  
  # Save as PDF
  ggsave(output_png, combined_plot, width = 16, height = 9)
}

# -------------------------
# COMMAND LINE INTERFACE
# -------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript admixture_plot.R <sample_order_file> <metadata_file> [output_png]\n", call. = FALSE)
}

sample_order_file <- args[1]
metadata_file <- args[2]
output_png <- ifelse(length(args) >= 3, args[3], "admixture_plots.png")

# Find all .Q files in the *current directory*
file_list <- list.files(pattern = "\\.Q$", full.names = TRUE)

if (length(file_list) == 0) {
  stop("No .Q files found in the current directory.\n", call. = FALSE)
}

# Run function
plot_admixture(file_list, sample_order_file, metadata_file, output_png)

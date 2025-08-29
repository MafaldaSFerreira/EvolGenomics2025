#!/usr/bin/env Rscript

# -------------------------
# LIBRARIES
# -------------------------
suppressPackageStartupMessages({
  library(tidyverse)
})

# -------------------------
# COMMAND LINE ARGUMENTS
# -------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript pca_plot.R <eigenvec_file> <eigenval_file> <metadata_file>\n", call. = FALSE)
}

eigenvec_file <- args[1]
eigenval_file <- args[2]
metadata_file <- args[3]

# -------------------------
# INPUTS
# -------------------------
metadata <- read.table(metadata_file, header = TRUE)

# -------------------------
# READ INPUT FILES
# -------------------------
pca <- read.table(eigenvec_file, header = FALSE)
eigenval <- scan(eigenval_file)

# Drop nuisance first column, set names to match
pca <- pca[, -1]
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca) - 1))

# Calculate % variance explained
pve <- eigenval / sum(eigenval) * 100

# Merge with metadata
pca_with_info <- merge(pca, metadata, by.x = "ind", by.y = "sample_id")

# -------------------------
# PLOT
# -------------------------
p <- ggplot(pca_with_info, aes(PC1, PC2, col = population, shape = species)) +
  geom_point(size = 2, alpha = 0.9)

p_maf5 <- p + labs(
  x = paste0("PC1 (", signif(pve[1], 3), "%)"),
  y = paste0("PC2 (", signif(pve[2], 3), "%)")
)

# -------------------------
# SAVE FIGURE
# -------------------------
ggsave(p_maf5, filename = "pca_plot.png", width = 6, height = 5)
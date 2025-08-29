#!/usr/bin/env Rscript

# ----------------------------
# Heterozygosity plots script
# ----------------------------

suppressPackageStartupMessages({
  library(tidyverse)
})

# COMMAND LINE ARGUMENTS ####
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: ./plot_heterozygosity.R <het_file> <metadata_file> <output_png>\n")
}

het_file     <- args[1]
metadata_file <- args[2]
output_png    <- args[3]

# READ DATA ####
het_tbl <- read.table(het_file, header = TRUE)
metadata <- read.table(metadata_file, header = TRUE)

# MERGE DATA ####
het_tbl_metadata <- merge(het_tbl, metadata,
                          by.x = "INDV", by.y = "sample_id")

# PLOTTING ####
pdf(output_png, width = 8, height = 6)

# Plot 1: Boxplot by population
print(
  ggplot(het_tbl_metadata) +
    geom_boxplot(aes(y = F, x = population_code)) +
    labs(title = "Inbreeding by Population",
         y = "F statistic",
         x = "Population")
)

# Plot 2: Proportion of observed homozygotes per population
print(
  ggplot(het_tbl_metadata) +
    geom_boxplot(aes(y = O.HOM./N_SITES, x = population_code)) +
    labs(title = "Observed Homozygote Proportion per Population",
         y = "O.HOM./N_SITES",
         x = "Population")
)

dev.off()

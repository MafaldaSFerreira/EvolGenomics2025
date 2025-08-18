# 

# LIBRARIES ####
library(tidyverse)

# INPUTS ####
# 1. Define your working directory:
setwd("~/Documents/Ferreira_SU/Repositories/EvolGenomics2025/pca/")

#-------------------------
# Sample info (used for all plots)
#-------------------------
metadata <- read.table("metadata/sample_metadata.txt", header=T)

#-------------------------
# File paths (your actual files)
#-------------------------
# Describe what these files contain:
eigenvec_file <- "02_results/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf3.var.10kbSNPs.chr1.27indv.eigenvec"
eigenval_file <- "02_results/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf3.var.10kbSNPs.chr1.27indv.eigenval"

# maf 5
eigenvec_file <- "02_results/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf5.var.10kbSNPs.chr1.27indv.eigenvec"
eigenval_file <- "02_results/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf5.var.10kbSNPs.chr1.27indv.eigenval"

#-------------------------
# Read input files
#-------------------------

pca <- read.table(eigenvec_file, header=F)
eigenval <- scan(eigenval_file)

# Drop nuisance first column, set names to match your code
pca <- pca[ , -1]
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# % variance explained
pve <- eigenval / sum(eigenval) * 100

# Merge with metadata + ensure stable factor order for colors
pca_with_info <- merge(pca, metadata, by.x = "ind", by.y = "sample_id")

# Plot

p <- ggplot(pca_with_info, aes(PC1, PC2, col = population, shape = species)) +
  geom_point(size = 2, alpha = 0.9)

# Let's add the information about the pve to each axis
p_maf3 <- p + labs(
    x = paste0("PC1 (", signif(pve[1], 3), "%)"),
    y = paste0("PC2 (", signif(pve[2], 3), "%)"))

p_maf5 <- p + labs(
  x = paste0("PC1 (", signif(pve[1], 3), "%)"),
  y = paste0("PC2 (", signif(pve[2], 3), "%)"))

gridExtra::grid.arrange(p_maf3, p_maf5)

ggsave(p_maf3, filename="figures/pca_maf3.pdf")
ggsave(p_maf5, filename="figures/pca_maf5.pdf")



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

if (length(args) < 1) {
  stop("Usage: Rscript cv_plot.R <cross_validation_file> [output_file]\n", call. = FALSE)
}

input_file <- args[1]
output_file <- ifelse(length(args) >= 2, args[2], "cross_validation.png")

# -------------------------
# READ INPUT FILE
# -------------------------
cross_validation_errors <- read.table(input_file, header = FALSE)
names(cross_validation_errors) <- c("K", "Error")

# -------------------------
# PLOT
# -------------------------
p <- ggplot(cross_validation_errors, aes(x = K, y = Error)) +
  geom_point() +
  geom_line() +
  labs(x = "K", y = "Cross-validation error")

# -------------------------
# SAVE FIGURE
# -------------------------
ggsave(p, filename = output_file, width = 6, height = 5)
# LIBRARIES ####
library(tidyverse)

# INPUTS ####
#Define your working directory:
setwd("~/Documents/Ferreira_SU/Repositories/EvolGenomics2025/02_population_structure_analysis/04_ADMIXTURE/")

#-------------------------
# Cross-validation Plot
#-------------------------
cross_validation_errors <- read.table("results/cross_validation_values.txt")
names(cross_validation_errors) <- c("K", "Error")

cross_validation_errors %>%
  ggplot() +
  geom_point(aes(x=K, y=Error))+
  geom_line(aes(x=K, y=Error))

#-------------------------
# Sample info (used for all plots)
#-------------------------

# A file with information about each individual, including the population and alternative names
metadata <- read.table("metadata/sample_metadata.txt", header=T)

# The order of the samples in the vcf and bed files.
sample_order <- read.table("metadata/sample_order.txt")

#-------------------------
# Admixture plots K=2
#-------------------------

# I will give you an example of how to plot K=2. Then you can use the function below to plot all K plots from a list of input files.

# We read the input file name
K2_file <- "results/chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.2.Q"
K2_tbl <- read.table(K2_file)

# Add a column with the names of the individuals
K2_tbl$indv<-sample_order$V1

# Add metadata to the table
K2_tbl_with_metadata <- merge(K2_tbl, metadata, by.x="indv", by.y="sample_id")

# Reorder the results by population for plotting. 
K2_tbl_with_metadata <- K2_tbl_with_metadata[order(K2_tbl_with_metadata$population_code),]

# Lock the order of the individuals for plotting.
K2_tbl_with_metadata$indv <- factor(K2_tbl_with_metadata$indv, levels=K2_tbl_with_metadata$indv)

# We need to change the format of the table for easier plotting with ggplots
K2_tbl_with_metadata_long <- K2_tbl_with_metadata %>% pivot_longer(cols=c("V1","V2"), names_to="K", values_to="admix_proportion")

# Plotting with ggplot:
ggplot(K2_tbl_with_metadata_long, aes(x=indv, y=admix_proportion, fill=K)) +
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        axis.title.x=element_blank())


library(cowplot)

plot_admixture <- function(files, sample_order, metadata, output_pdf="admixture_plots_grid.pdf") {
  
  plots <- list()
  
  for (f in files) {
    # Read Q file
    Q_tbl <- read.table(f)
    
    # Add individuals
    Q_tbl$indv <- sample_order$V1
    
    # Merge with metadata
    Q_tbl_with_metadata <- merge(Q_tbl, metadata, by.x="indv", by.y="sample_id")
    
    # Reorder individuals by population
    Q_tbl_with_metadata <- Q_tbl_with_metadata[order(Q_tbl_with_metadata$population_code), ]
    Q_tbl_with_metadata$indv <- factor(Q_tbl_with_metadata$indv, 
                                       levels=Q_tbl_with_metadata$indv)
    
    # Convert to long format
    Q_tbl_long <- Q_tbl_with_metadata %>%
      pivot_longer(cols = starts_with("V"), 
                   names_to = "K", 
                   values_to = "admix_proportion")
    
    # Make plot
    p <- ggplot(Q_tbl_long, aes(x=indv, y=admix_proportion, fill=K)) +
      geom_bar(position="stack", stat="identity") +
      theme_classic() +
      theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
            axis.title.x=element_blank()) +
      ggtitle(basename(f))
    
    plots[[f]] <- p
  }
  
  # Combine all plots in a grid
  combined_plot <- plot_grid(plotlist = plots, ncol = 2)  # change ncol if you want more/less per row
  
  # Save as PDF
  ggsave(output_pdf, combined_plot, width = 16, height = 9)
}


file_list <- list.files("results", pattern="*.Q", full.names=TRUE)

plot_admixture(file_list, sample_order, metadata, "admixture_plots.pdf")      

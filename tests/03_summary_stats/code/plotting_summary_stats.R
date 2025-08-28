# LIBRARIES ####
library(tidyverse)

# INPUTS ####
# 1. Define your working directory:
setwd("~/Documents/Ferreira_SU/Repositories/EvolGenomics2025/summary_stats/")

# wFST and meanFst

wFst <- read.table("results/wFst.out", header = F)
colnames(wFst) <- c("contrast", "wFst")

meanFst <- read.table("results/meanFst.out", header = F)
colnames(meanFst) <- c("contrast", "meanFst")

# Let's combine both values into one single table
Fst_tbl <- merge(meanFst, wFst, by="contrast")

# Let's create two columns containing the names of pop1 and pop2
Fst_tbl$pop1 <- str_split_fixed(Fst_tbl$contrast, "_vs_", 2)[,1]
Fst_tbl$pop2 <- str_split_fixed(Fst_tbl$contrast, "_vs_", 2)[,2]


# let's plot the Fst values as a heatmap
ggplot(Fst_tbl) +
  geom_tile(aes(x=pop1, y=pop2, fill=wFst))+
  geom_text(aes(x=pop1, y=pop2, label=wFst), color="white")
  

ggplot(Fst_tbl) +
  geom_tile(aes(x=pop1, y=pop2, fill=meanFst))+
  geom_text(aes(x=pop1, y=pop2, label=meanFst), color="white")


# Plotting heterozygosity:

het_tbl <- read.table("results/heterozygosity.out.het", header=T)

# let's plot by individual:

ggplot(het_tbl) +
  geom_boxplot(aes(y=F))

# here, it will be nice to have metadata information added to the table
# so that we can plot the results by population.
metadata <- read.table("metadata/sample_metadata.txt", header=T)

het_tbl_metadata <- merge(het_tbl, metadata, by.x="INDV", by.y="sample_id")

ggplot(het_tbl_metadata) +
  geom_boxplot(aes(y=F, x=population_code))

# We can also plot the proportion of observed homozygote sites 
# per population. What do you observe? How does that relate to the values above?
ggplot(het_tbl_metadata) +
  geom_boxplot(aes(y=O.HOM./N_SITES, x=population_code))




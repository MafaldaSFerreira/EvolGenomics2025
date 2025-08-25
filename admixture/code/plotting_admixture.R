# LIBRARIES ####
library(tidyverse)

# INPUTS ####
# 1. Define your working directory:
setwd("~/Documents/Ferreira_SU/Repositories/EvolGenomics2025/admixture/")

#-------------------------
# Sample info (used for all plots)
#-------------------------
metadata <- read.table("metadata/sample_metadata.txt", header=T)
sample_order <- read.table("metadata/sample_order.txt")

#-------------------------
# Cross-validation Plot
#-------------------------
cross_validation_errors <- read.table("results/cross_validation_to_plot.txt")
names(cross_validation_errors) <- c("K", "Error")

cross_validation_errors %>% #filter(K<=5) %>%
ggplot() +
  geom_point(aes(x=K, y=Error))+
  geom_line(aes(x=K, y=Error))

# The best K seems to be 2 or 3. So let's plot those:



#-------------------------
# Admixture plots K=2
#-------------------------
K2_file <- "results/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf5.var.10kbSNPs.chr1.27indv.2.Q"
K2_tbl <- read.table(K2_file)

K2_tbl$indv<-sample_order$V1

K2_tbl_with_metadata <- merge(K2_tbl, metadata, by.x="indv", by.y="sample_id")

# Reorder by population
K2_tbl_with_metadata <- K2_tbl_with_metadata[order(K2_tbl_with_metadata$population_code),]

K2_tbl_with_metadata$indv <- factor(K2_tbl_with_metadata$indv, levels=K2_tbl_with_metadata$indv)

K2_tbl_with_metadata %>%
  pivot_longer(cols=c("V1","V2"), 
               names_to="K", values_to="admix_proportion") %>%
  ggplot(aes(x=indv, y=admix_proportion, fill=K)) +
  geom_bar(position="stack", stat="identity") +
    theme_classic() +
    theme(axis.title.y=element_text(size=6),
          axis.text.y=element_text(size=6),
          axis.text.x=element_text(size=6,angle=90, vjust=0.5, hjust=1),
          axis.title.x=element_blank())


#-------------------------
# Admixture plots K=3
#-------------------------
K3_file <- "results/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf5.var.10kbSNPs.chr1.27indv.3.Q"
K3_tbl <- read.table(K3_file)

K3_tbl$indv<-sample_order$V1

K3_tbl_with_metadata <- merge(K3_tbl, metadata, by.x="indv", by.y="sample_id")

# Reorder by population
K3_tbl_with_metadata <- K3_tbl_with_metadata[order(K3_tbl_with_metadata$population_code),]

K3_tbl_with_metadata$indv <- factor(K3_tbl_with_metadata$indv, levels=K3_tbl_with_metadata$indv)

K3_tbl_with_metadata %>%
  pivot_longer(cols=c("V1","V2", "V3"), 
               names_to="K", values_to="admix_proportion") %>%
  ggplot(aes(x=indv, y=admix_proportion, fill=K)) +
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  theme(axis.title.y=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6,angle=90, vjust=0.5, hjust=1),
        axis.title.x=element_blank())

#-------------------------
# Admixture plots K=4
#-------------------------
K4_file <- "results/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf5.var.10kbSNPs.chr1.27indv.4.Q"
K4_tbl <- read.table(K4_file)

K4_tbl$indv<-sample_order$V1

K4_tbl_with_metadata <- merge(K4_tbl, metadata, by.x="indv", by.y="sample_id")

# Reorder by population
K4_tbl_with_metadata <- K4_tbl_with_metadata[order(K4_tbl_with_metadata$population_code),]

K4_tbl_with_metadata$indv <- factor(K4_tbl_with_metadata$indv, levels=K4_tbl_with_metadata$indv)

K4_tbl_with_metadata %>%
  pivot_longer(cols=c("V1","V2", "V3", "V4"), 
               names_to="K", values_to="admix_proportion") %>%
  ggplot(aes(x=indv, y=admix_proportion, fill=K)) +
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  theme(axis.title.y=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6,angle=90, vjust=0.5, hjust=1),
        axis.title.x=element_blank())

#-------------------------
# Admixture plots K=5
#-------------------------
K5_file <- "results/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf5.var.10kbSNPs.chr1.27indv.5.Q"
K5_tbl <- read.table(K5_file)

K5_tbl$indv<-sample_order$V1

K5_tbl_with_metadata <- merge(K5_tbl, metadata, by.x="indv", by.y="sample_id")

# Reorder by population
K5_tbl_with_metadata <- K5_tbl_with_metadata[order(K5_tbl_with_metadata$population_code),]

K5_tbl_with_metadata$indv <- factor(K5_tbl_with_metadata$indv, levels=K5_tbl_with_metadata$indv)

K5_tbl_with_metadata %>%
  pivot_longer(cols=c("V1","V2", "V3", "V4", "V5"), 
               names_to="K", values_to="admix_proportion") %>%
  ggplot(aes(x=indv, y=admix_proportion, fill=K)) +
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  theme(axis.title.y=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6,angle=90, vjust=0.5, hjust=1),
        axis.title.x=element_blank())

# Pophelper

# -------------------------
# 1. Collect all Q files
# -------------------------
# Adjust the pattern to match your Q files. 
# Example: run1/K2.Q, run2/K2.Q, etc.

# Read CV errors
cv <- read.table("results/cv_errors.txt", header = FALSE)
colnames(cv) <- c("K", "CVerror")

# Summarize mean & sd per K
cv_summary <- cv %>%
  group_by(K) %>%
  summarise(
    mean_cv = mean(CVerror),
    sd_cv   = sd(CVerror),
    .groups = "drop"
  )

print(cv_summary)

# Plot
p <- ggplot(cv_summary, aes(x = K, y = mean_cv)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_cv - sd_cv, ymax = mean_cv + sd_cv), width = 0.2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "ADMIXTURE cross-validation error",
    x = "Number of clusters (K)",
    y = "Mean CV error ± SD"
  )


qfiles <- list.files(path="admixture_results", pattern="\\.Q$", recursive=TRUE, full.names=TRUE)
names(qfiles) <- gsub("/","",qfiles)

#Read them in
slist <- readQ(qfiles)

names(slist) <- gsub("/","",qfiles)

# -------------------------
# 2. Summarize runs for each K
# -------------------------
# This will cluster replicate runs at each K and produce 
# an alignment of clusters (so cluster 1 in run1 ~ cluster 1 in run2, etc.)
# You can increase rep=10 if you want more iterations of alignment
aligned <- alignK(slist )


plotQ(aligned)

pdf("admixture_barplots.pdf", width=10, height=5)
for(k in sort(unique(sapply(slist, function(x) ncol(x))))){
  plist <- aligned[[paste0("K",k)]]
  if(!is.null(plist)){
    plotQ(plist,
          showindlab=FALSE,
          sortind="all",   # sorts individuals within each cluster
          exportpath=NULL, # don’t auto-save PNGs, we’re inside pdf()
          imgtype="pdf")
  }
}
dev.off()


# This is too complicated, I think I will have a first part where they run jsut 1 replicate
# to understood how to do it, and then I can give them further code to make this work for several 
#runs? 



  
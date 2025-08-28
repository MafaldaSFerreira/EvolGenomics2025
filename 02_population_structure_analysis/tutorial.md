# Population Structure Analysis 

In this tutorial, we will cover how to use a population vcf file to study population structure in our dataset using Principal Component Analysis, Admixture and summary statistics (Fst)


## Principal Component Analysis with PLINK



We will use `plink` to perform a PCA analysis with our set of biallelic, independent SNPs. We will place the output files in a separate directory, to keep our working directory organized.

To generate the pca, we run plink with the flag `--pca`. 

Note that when specifying the input files with flag `--bfile`, you only specify the prefix that is common to `.bed`, `.bim` and `.fam` files.

run_plink_pca.sh
~~~bash
#!/bin/bash
#SBATCH -A naiss2025-22-172
#SBATCH --job-name=plink
#SBATCH --output=logs/plink_%A_%a.out
#SBATCH --error=logs/plink_%A_%a.err
#SBATCH -p shared
#SBATCH -n 6
#SBATCH -N 1
#SBATCH --time=1-00:00:00 
#SBATCH --mail-type=all
#SBATCH --mail-user=sferreira.mafalda@gmail.com

# load necessary software
ml load bioinfo-tools plink

# Define the path to your input file:
INPUT_DIR="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data/02_filtered_vcfs_for_analysis"

# Define the path where you want to store the output file:
OUTPUT_DIR="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data/03_PCA"

# This command with create the output directory, if it doesn't exist:
mkdir -p ${OUTPUT_DIR}

# Run plink
plink --bfile ${INPUT_DIR}/chr1.27indvs.QUALFilters.POPFilters.10kbSNPs --pca --out ${OUTPUT_DIR}/chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.PCA
~~~

In the output folder, you will find four files.

~~~
# A file containing the percentage variation explained by the first 20 PC axis.
chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.PCA.eigenval

# A matrix file containing per individual  PC loadings for PC1 to PC20. 
chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.PCA.eigenvec  

# List of individuals
chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.PCA.nosex

# Log of the run
chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.PCA.log
~~~

We want to use `.eigenvec` and `.eigenval` to plot in R. So download them to your local folder.

~~~R
#Plotting files in R

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

# Calculate the % variance explained
pve <- eigenval / sum(eigenval) * 100

# Merge with metadata + ensure stable factor order for colors
pca_with_info <- merge(pca, metadata, by.x = "ind", by.y = "sample_id")

# Plot
p <- ggplot(pca_with_info, aes(PC1, PC2, col = population, shape = species)) +
  geom_point(size = 2, alpha = 0.9)

# Let's add the information about the pve to each axis
p_maf5 <- p + labs(
  x = paste0("PC1 (", signif(pve[1], 3), "%)"),
  y = paste0("PC2 (", signif(pve[2], 3), "%)"))

ggsave(p_maf5, filename="figures/pca.pdf")
~~~

>Q. How do populations separate in the Principal Component Analysis? What does the distance between suggest about their differentiation?


## Estimating Ancestry using ADMIXTURE

In complement to PCA, we will use the software `ADMIXTURE` to run an ancestry analysis. 

We will run the software 6 times, each time changing the expected value of K , or ancestral populations, from 1 to 6. The reason to choose 6 comes from our a priori knowledge of the dataset. We know we have sampled individuals from 4 regions so that might be the most likely number of ancestral populations, but we can run admixture until K=6 to understand if any individual show signs of admixed ancestry. 

To simplify our task and plotting, and because our SNP dataset is very simple, we are going to run admixture only for one iteration per K. However, it is common practice to run several iterations of admixture (10 to 100s). 

~~~
    # Run admixture with 10 fold cross-validation 
    --cv 10

    # We pipe the log output to a separate file, so that we can access the likelihood and cross validation error for each run
    tee log_k${k}.out
~~~


~~~bash
#!/bin/bash
#SBATCH -A naiss2025-22-172
#SBATCH --job-name=admixture
#SBATCH --output=logs/admixture_%A_%a.out
#SBATCH --error=logs/admixture_%A_%a.err
#SBATCH -p shared
#SBATCH -n 6
#SBATCH -N 1
#SBATCH --time=1-00:00:00 
#SBATCH --mail-type=all
#SBATCH --mail-user=sferreira.mafalda@gmail.com

# load necessary software
ml load bioinfo-tools ADMIXTURE

# Define the path to your input file:
INPUT_DIR="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data/02_filtered_vcfs_for_analysis"

# Define the path where you want to store the output file:
OUTPUT_DIR="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data/04_ADMIXTURE"

# This command with create the output directory, if it doesn't exist:
mkdir -p ${OUTPUT_DIR}

# Admixture will generate the output results in the current directory, so move to the output directory:
cd ${OUTPUT_DIR}

# Let's run admixture for 1 to 6 values of K:
for k in $(seq 1 6); 
    do admixture --cv=10 ${INPUT_DIR}/chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.bed ${k} | tee log_k${k}.out; 
done;
~~~

After ADMIXTURE has run, navigate inside the OUTPUT_DIR file, where all the outputs are stored. Here you will find three files per run:

~~~
# P (the allele frequencies of the inferred ancestral populations)
chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.1.P

# Q (the ancestry fractions); which we use for plotting:
chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.1.Q

# Log output of the run, which contains the likelihood and cross validation error.
log_k1.out
~~~

Before we plot the admixture results, we will compile the cross validation error results to estimate the best value of K in our dataset. We can run the code below to perform this task.

~~~bash
cd 04_ADMIXTURE/
grep -h CV log*.out 

#CV error (K=1): 0.81212
#CV error (K=2): 0.56319
#CV error (K=3): 0.58352
#CV error (K=4): 0.61644
#CV error (K=5): 0.75005
#CV error (K=6): 0.80546

grep -h CV log*.out | cut -d" " -f3,4 | sed 's/(K=//g' | sed 's/)://g' | tr " " "\t" > cross_validation_values.txt
~~~

Use `cat` to make sure the file `cross_validation_values.txt` is well formatted.

Switch to R_studio and plot the `cross_validation_values.txt`.

~~~R
# LIBRARIES ####
library(tidyverse)

# INPUTS ####
#Define your working directory:
setwd("~/Documents/Ferreira_SU/Repositories/EvolGenomics2025/admixture/")

#-------------------------
# Cross-validation Plot
#-------------------------
cross_validation_errors <- read.table("results/cross_validation_to_plot.txt")
names(cross_validation_errors) <- c("K", "Error")

cross_validation_errors %>%
ggplot() +
  geom_point(aes(x=K, y=Error))+
  geom_line(aes(x=K, y=Error))

~~~

> Q. What seems to be the best K value? Hint: A good value of K will exhibit a low cross-validation error compared to other K values. Hint2: Also take into consideration the PCA results and you're knowledge of the dataset.


~~~R
#(continues from code above)
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
K2_tbl_with_metadata_long <- pivot_longer(cols=c("V1","V2"), names_to="K", values_to="admix_proportion")

# Plotting with ggplot:
ggplot(K2_tbl_with_metadata_long, aes(x=indv, y=admix_proportion, fill=K)) +
  geom_bar(position="stack", stat="identity") +
    theme_classic() +
    theme(axis.title.y=element_text(size=6),
          axis.text.y=element_text(size=6),
          axis.text.x=element_text(size=6,angle=90, vjust=0.5, hjust=1),
          axis.title.x=element_blank())
~~~

We could save the plot we just generated in the following way:

~~~R

K2_plot <- ggplot(K2_tbl_with_metadata_long, aes(x=indv, y=admix_proportion, fill=K)) +
  geom_bar(position="stack", stat="identity") +
    theme_classic() +
    theme(axis.title.y=element_text(size=6),
          axis.text.y=element_text(size=6),
          axis.text.x=element_text(size=6,angle=90, vjust=0.5, hjust=1),
          axis.title.x=element_blank())

ggsave(K2_plot, filename="figures/K2_plot.pdf", width=10, height=5)
~~~

Great, we learned how to do one plot. 

I have generated a function below that allows us to plot all Ks at once and outputs a pdf file.

~~~R
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

# Read the list of .Q files.
file_list <- list.files("results", pattern="*.Q", full.names=TRUE)

# Generate one plot per K and output a single pdf page with the plots.
plot_admixture(file_list, sample_order, metadata, "admixture_plots.pdf")      
    
~~~

## Genetic differentiation between populations and levels of diveristy 

The previous analysis alowed us to confirm the number of populations we have in our dataset. We learned that we likely have 4 populations, with a low level of admixture between then. With this information, we will now calculate genetic differentiation (Fst) between populations to understand how differentiated they are. We will also assess levels of genetic diveristy (heterozygosity) and inbreeding (FIS) to learn more about the demography of the populations. 

We will again use `vcftools` to make the calculations we need. To calculate Fst among populations, we will use the same set of independent SNPs we used for PCA and ADMIXTURE analysis. We will perform calculations in a pairwise fassion between populations. We will several input files contained in `metadata/` which contain a list of individual names per population.

~~~bash
ls metadata/
FRA.txt  GER.txt  HEL.txt  IND.txt
~~~

~~~bash
#!/bin/bash
#SBATCH -A naiss2025-22-172
#SBATCH --job-name=FST
#SBATCH --output=logs/FST_%A_%a.out
#SBATCH --error=logs/FST_%A_%a.err
#SBATCH -p shared
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --time=01:00:00 
#SBATCH --mail-type=all
#SBATCH --mail-user=sferreira.mafalda@gmail.com

# Load the necessary software
ml load bioinfo-tools vcftools

# Define the path to your input file:
INPUT_DIR="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data/02_filtered_vcfs_for_analysis"

# Metadata directory
METADATA_DIR="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data/metadata"

# Define the path where you want to store the output file:
OUTPUT_DIR="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data/05_FST"

# This command with create the output directory, if it doesn't exist:
mkdir -p ${OUTPUT_DIR}

cd ${OUTPUT_DIR}

# Run vcftools for each population pair
vcftools --gzvcf ${INPUT_DIR}/chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.vcf.gz --weir-fst-pop ${METADATA_DIR}/FRA.txt --weir-fst-pop ${METADATA_DIR}/IND.txt --out FRA_vs_IND

vcftools --gzvcf ${INPUT_DIR}/chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.vcf.gz --weir-fst-pop ${METADATA_DIR}/FRA.txt --weir-fst-pop ${METADATA_DIR}/GER.txt --out FRA_vs_GER

vcftools --gzvcf ${INPUT_DIR}/chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.vcf.gz --weir-fst-pop ${METADATA_DIR}/FRA.txt --weir-fst-pop ${METADATA_DIR}/HEL.txt --out FRA_vs_HEL

vcftools --gzvcf ${INPUT_DIR}/chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.vcf.gz --weir-fst-pop ${METADATA_DIR}/IND.txt --weir-fst-pop ${METADATA_DIR}/GER.txt --out IND_vs_GER

vcftools --gzvcf ${INPUT_DIR}/chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.vcf.gz --weir-fst-pop ${METADATA_DIR}/IND.txt --weir-fst-pop ${METADATA_DIR}/HEL.txt --out IND_vs_HEL

vcftools --gzvcf ${INPUT_DIR}/chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.vcf.gz --weir-fst-pop ${METADATA_DIR}/GER.txt --weir-fst-pop ${METADATA_DIR}/HEL.txt --out GER_vs_HEL
~~~

After the commands has run, you will find a set of `.log` files which will contain the Fst calculations per population.

Use `grep` to extract them all comparisons at the same time.

~~~
grep "Weir and Cockerham mean Fst estimate:" *log | cut -f1,3 -d":" | sed 's/.log:/\t/' > meanFst.out
grep "Weir and Cockerham weighted Fst estimate:" *log | cut -f1,3 -d":" | sed 's/.log:/\t/' > wFst.out
~~~

Use `cat` to inspect the file.

~~~bash
#!/bin/bash
#SBATCH -A naiss2025-22-172
#SBATCH --job-name=het
#SBATCH --output=logs/het_%A_%a.out
#SBATCH --error=logs/het_%A_%a.err
#SBATCH -p shared
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --time=01:00:00 
#SBATCH --mail-type=all
#SBATCH --mail-user=sferreira.mafalda@gmail.com

# Load the necessary software
ml load bioinfo-tools vcftools

# Define the path to your input file:
INPUT_DIR="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data/02_filtered_vcfs_for_analysis"

# Define the path where you want to store the output file:
OUTPUT_DIR="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data/05_FST"

# This command with create the output directory, if it doesn't exist:
mkdir -p ${OUTPUT_DIR}

cd ${OUTPUT_DIR}

vcftools --gzvcf ${INPUT_DIR}/chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.vcf.gz --het --out het_fis_calculations
~~~

Here, the output file we are interested is `het_fis_calculations.het`.


~~~R
# LIBRARIES ####
library(tidyverse)

# INPUTS ####
# 1. Define your working directory:
setwd("~/Documents/Ferreira_SU/Repositories/EvolGenomics2025/summary_stats/")

# wFST and meanFst
wFst <- read.table("results/wFst.out", header = F)
meanFst <- read.table("results/meanFst.out", header = F)

# Add column names to the tables
colnames(wFst) <- c("contrast", "wFst")
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
~~~




~~~R
# LIBRARIES ####
library(tidyverse)

# INPUTS ####
# 1. Define your working directory:

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
~~~
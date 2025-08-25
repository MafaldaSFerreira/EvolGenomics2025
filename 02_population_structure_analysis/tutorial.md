# Population Structure Analysis 

In this tutorial, we will cover how to use a population vcf file to study population structure in our dataset using Principal Component Analysis, Admixture and summary statistics (Fst)

## Applying Population filters to SNP dataset
We need allele frequency data for these analyses, which we calculate from SNP data. To ensure condifence in our calculations, we will keep only sites with two alleles (i.e., biallelic sites) and for which the minor allele frequency is higher than 0.05. This ensures we are not calculating allele frequencies from rare variants that could be technical errors generated during sequencing. We also only want to keep sites where we have information for a minimum number of individuals.

The tool we are using today to perform this filtering is `vcftools`. 

Let's breakdown the filters in the vcftools command below:
~~~
    # Exclude sites where more than 20% of the individuals have missing genotypes
    --max-missing 0.2 \ 
    # Exclude sites with less than 2 alleles
    --min-alleles 2 \ 
    # Exclude sites with more than 2 alleles
    --max-alleles 2 \ 
    # Exclude sites where the minor allele has a frequency lower than 0.05
    --maf 0.05 \ 
    # Output a vcf file
    --recode --recode-INFO-all --stdout | bgzip -c > ${OUTPUT_DIR}/chr1.27indvs.QUALFilters.POPFilters.vcf.gz
~~~

Now we can run this as a bash script on the cluster:
~~~bash
#!/bin/bash
#SBATCH -A naiss2025-22-172
#SBATCH --job-name=vcftools
#SBATCH --output=logs/vcftools_%A_%a.out
#SBATCH --error=logs/vcftools_%A_%a.err
#SBATCH -p shared
#SBATCH -n 6
#SBATCH -N 1
#SBATCH --time=1-00:00:00 
#SBATCH --mail-type=all
#SBATCH --mail-user=sferreira.mafalda@gmail.com

# load necessary software
ml load bioinfo-tools vcftools

# Define the path to your input vcf file:
INPUT_DIR="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data/01_filtered_vcfs"

# Define the path where you want to store the output file:
OUTPUT_DIR="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data/02_filtered_vcfs_for_analysis"

# Run vcftools_
vcftools --gzvcf ${INPUT_DIR}/chr1.27indvs.QUALFilters.vcf.gz --max-missing 0.2 --min-alleles 2 --max-alleles 2 --maf 0.05 --recode --recode-INFO-all --stdout | bgzip -c > ${OUTPUT_DIR}/chr1.27indvs.QUALFilters.POPFilters.vcf.gz
~~~


Another important filter we will do for the population structure analysis is one where we will thin our SNP dataset so that we keep SNPs that are least 10kb from each other. This allows us to keep variants that may be evolving approximately independetly from each other, since they will be in separate recombination blocks. SNPs that are linked to each other will give the same information about neutral population demography. We will create a separate vcf file here, since we might still want to do analysis with all SNPs later on, like calculating heterozygosity.

~~~
# Keep SNPs that are at least 10kb from each other
    --thin 10000
~~~


~~~bash
#!/bin/bash
#SBATCH -A naiss2025-22-172
#SBATCH --job-name=vcftools
#SBATCH --output=logs/thin_%A_%a.out
#SBATCH --error=logs/thin_%A_%a.err
#SBATCH -p shared
#SBATCH -n 6
#SBATCH -N 1
#SBATCH --time=1-00:00:00 
#SBATCH --mail-type=all
#SBATCH --mail-user=sferreira.mafalda@gmail.com

# load necessary software
ml load bioinfo-tools vcftools

# Define the path to your input vcf file:
INPUT_DIR="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data/02_filtered_vcfs_for_analysis"

# Define the path where you want to store the output file:
OUTPUT_DIR="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data/02_filtered_vcfs_for_analysis"

vcftools --gzvcf ${INPUT_DIR}/chr1.27indvs.QUALFilters.POPFilters.vcf.gz --thin 10000 --recode --recode-INFO-all --stdout | bgzip -c > ${OUTPUT_DIR}/chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.vcf.gz  
~~~

Finally, we will need plink input files which are in .bed .bim and .fam formats. Since we will need them both for PCA and admixture analysis, let's generate this plink input files in the same 02_filtered_vcfs_for_analysis folder.

~~~
    # causes both family and within-family IDs to be set to the sample ID
    --double-id
    # to allow non human chromosome codes
    --allow-extra-chr
    # adds IDs to variants
    --set-missing-var-ids @:#
    # make the bed input files
    --make-bed
    # name the bed files
    --out ${OUTPUT_DIR}/chr1.27indvs.QUALFilters.POPFilters.10kbSNPs
~~~

run_plink_makebed.sh
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
OUTPUT_DIR="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data/02_filtered_vcfs_for_analysis"

plink --vcf ${INPUT_DIR}/chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --out ${OUTPUT_DIR}/chr1.27indvs.QUALFilters.POPFilters.10kbSNPs

~~~

To facilitate our analysis, I have placed a set of metadata files on the server. 

~~~bash
ls metadata/
FRA.txt  GER.txt  HEL.txt  IND.txt  populations.txt  samples.txt
~~~

`samples.txt` : contains the list of all 27 individuals in the order they appear in the vcf file

`populations.txt` : contains the population assignment of each individual, in the order the individuals are in the vcf file

`FRA.txt`, `GER.txt`, `HEL.txt` and `IND.txt`: files containing the list of individuals from each population.

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


We will now move to R to plot the output files of the PCA:
~~~R
#Plotting files in R

~~~




## Estimating Ancestry using ADMIXTURE

We will use the software `ADMIXTURE` to run an ancestry analysis. 

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

After admixture has run, navigate inside the OUTPUT_DIR file, where all the outputs are stored.

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

ggsave(K2_plot, filename="K2_plot.pdf", width=10, height=5)
~~~

Great , we learned how to do one plot. 

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

> Q. What is your conclusion about the best number of K?
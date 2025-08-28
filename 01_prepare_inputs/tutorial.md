## Applying Quality filters to SNP dataset


We are working with a vcf file for chromosome 1, containing called genotypes for 27 individuals. It is good practice to apply filters to the called genotypes before continuing with formal analysis. We will have different levels of confidence for the accuracy of the genotypes that were called depending on depth of coverage (how many reads map to a particular position), genotype quality assigned by the software we used to call genotypes, if variants are SNPs or INDELs, etc.

So let's perform some general filtering. For this, we will use the program bcftools. Open your terminal and log-in to the server. The vcf file we are going to filter is the one you have generated already in the previous classes. 

You could take a PhD on how to filter genotypes. Here I apply best practices.

The first step applied across rows in the vcf file, includes:
- QUAL >= 30 : Keep genotypes with genotype quality above 30
- MQ >= 30: The mapping quality of the reads mapping to this position is above 30.
- INFO/DP < 1900: The average depth across all individuals (the INFO field) is less than 1900. This prevents us from including calls in repetitive regions of the genome.

Second step, applied to each genotype in a single row, includes:
- TYPE="snp": Keep only SNPs
- FMT/QG>=20 : The individual's genotype quality is >= 20
- FMT/DP>5 : Depth for the individual's genotype call is above 5.

Final filter, excludes:
- ALT=* : Exclude sites where the alternative allele is any (*) nucleotide
- TYPE~"indel": Exclude indels
- REF="N": Exclude sites where the reference allele is missing (or N)

~~~bash
#!/bin/bash
#SBATCH -A naiss2025-22-172
#SBATCH --job-name=filtering
#SBATCH --output=logs/filtering_%A_%a.out
#SBATCH --error=logs/filtering_%A_%a.err
#SBATCH -p shared
#SBATCH -n 6
#SBATCH -N 1
#SBATCH --time=1-00:00:00 
#SBATCH --mail-type=all
#SBATCH --mail-user=sferreira.mafalda@gmail.com

WD="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data"
INPUTWD=${WD}"/00_unfiltered_vcfs/subset_vcf"

ml load bioinfo-tools bcftools

cd $INPUTWD

# bcftools filter for quality, mapping quality and maximum depth
bcftools filter -i 'QUAL>=30 && MQ>=30 && INFO/DP<1900' -Oz -o chr1.27indv.QUAL30_MQ30_depthFilter1900.vcf.gz chr1.27indv.vcf.gz

# Sample specific filters:
bcftools filter -S . -i 'TYPE="snp" & FMT/GQ>=20 & FMT/DP>5' -Oz -o chr1.27indv.QUAL30_MQ30_depthFilter1900_GQ20_DP5.vcf.gz chr1.27indv.QUAL30_MQ30_depthFilter1900.vcf.gz

# Finally, select variants that have PASS the above filters, are not indels, don't have * ALT alleles and where the reference is not a missing genotype "N"
bcftools view -f PASS -e 'ALT="*" | TYPE~"indel" | ref="N"' -O z -o chr1.27indvs.QUALFilters.vcf.gz chr1.27indv.QUAL30_MQ30_depthFilter1900_GQ20_DP5.vcf.gz

bcftools index chr1.27indvs.QUALFilters.vcf.gz
bcftools stats -s - chr1.27indvs.QUALFilters.vcf.gz > chr1.27indvs.QUALFilters.stats
~~~

## Applying Population filters to SNP dataset
We need allele frequency data for our analayses, which we calculate from SNP data. To ensure condifence in our calculations, we will keep only sites with two alleles (i.e., biallelic sites) and for which the minor allele frequency is higher than 0.05. This ensures we are not calculating allele frequencies from rare variants that could be technical errors generated during sequencing. We also only want to keep sites where we have information for a minimum number of individuals.

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

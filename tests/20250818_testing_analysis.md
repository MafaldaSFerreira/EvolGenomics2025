# 2025-08-18

Preparing the contents for class

I want to test which regions of the genome might be good for plotting

We might have to later on change the way we filter the dataset

Let's subset the filtered vcf file 

run_bcftools_tests.sh
~~~bash
#!/bin/bash
#SBATCH -A naiss2025-22-172
#SBATCH --job-name=filtering
#SBATCH --output=filtering_%A_%a.out
#SBATCH --error=filtering_%A_%a.err
#SBATCH -p shared
#SBATCH -n 6
#SBATCH -N 1
#SBATCH --time=1-00:00:00 
#SBATCH --mail-type=all
#SBATCH --mail-user=sferreira.mafalda@gmail.com

ml load bioinfo-tools bcftools

WD="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data"
INPUTWD=${WD}"/01_filtered_vcfs"
METADWD=${WD}"/metadata"

cd ${INPUTWD}

bcftools stats -s - AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.vcf.gz > AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.vcf.stats

bcftools view -O z -o AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.chr1.vcf.gz AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.vcf.gz 1 

bcftools view -O z -o AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.chr7.vcf.gz AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.vcf.gz 7 

bcftools index AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.chr1.vcf.gz
bcftools index AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.chr7.vcf.gz


bcftools view --samples-file ${METADWD}/samples.txt -O z -o AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.chr7.27indv.vcf.gz AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.chr7.vcf.gz

bcftools view --samples-file ${METADWD}/samples.txt -O z -o AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.chr1.27indv.vcf.gz AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.chr1.vcf.gz
~~~


run_bcftools_tests.sh
~~~bash
#!/bin/bash
#SBATCH -A naiss2025-22-172
#SBATCH --job-name=filtering
#SBATCH --output=filtering_%A_%a.out
#SBATCH --error=filtering_%A_%a.err
#SBATCH -p shared
#SBATCH -n 6
#SBATCH -N 1
#SBATCH --time=1-00:00:00 
#SBATCH --mail-type=all
#SBATCH --mail-user=sferreira.mafalda@gmail.com

ml load bioinfo-tools bcftools

WD="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data"
INPUTWD=${WD}"/01_filtered_vcfs"
METADWD=${WD}"/metadata"

cd ${INPUTWD}

bcftools stats -s - AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.vcf.gz > AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.vcf.stats

bcftools view -O z -o AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.chr1.vcf.gz AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.vcf.gz 1 

bcftools view -O z -o AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.chr7.vcf.gz AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.vcf.gz 7 

bcftools index AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.chr1.vcf.gz
bcftools index AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.chr7.vcf.gz


bcftools view --samples-file ${METADWD}/samples.txt -O z -o AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.chr7.27indv.vcf.gz AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.chr7.vcf.gz

bcftools view --samples-file ${METADWD}/samples.txt -O z -o AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.chr1.27indv.vcf.gz AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.chr1.vcf.gz
~~~


run_bcftools_tests_unfiltered.sh
~~~bash
#!/bin/bash
#SBATCH -A naiss2025-22-172
#SBATCH --job-name=filtering
#SBATCH --output=filtering_%A_%a.out
#SBATCH --error=filtering_%A_%a.err
#SBATCH -p shared
#SBATCH -n 6
#SBATCH -N 1
#SBATCH --time=1-00:00:00 
#SBATCH --mail-type=all
#SBATCH --mail-user=sferreira.mafalda@gmail.com

ml load bioinfo-tools bcftools

WD="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data"
INPUTWD=${WD}"/00_unfiltered_vcfs"
METADWD=${WD}"/metadata"

cd ${INPUTWD}

bcftools view --samples-file ${METADWD}/samples.txt -O z -o chr1.27indv.vcf.gz 1.vcf.gz 
bcftools view --samples-file ${METADWD}/samples.txt -O z -o chr7.27indv.vcf.gz 7.vcf.gz 

bcftools index chr1.27indv.vcf.gz
bcftools index chr7.27indv.vcf.gz

bcftools stats -s - chr1.27indv.vcf.gz > chr1.27indv.vcf.stats
bcftools stats -s - chr7.27indv.vcf.gz > chr7.27indv.vcf.stats

~~~

run_bcftools_hard_filter.sh
~~~bash
#!/bin/bash
#SBATCH -A naiss2025-22-172
#SBATCH --job-name=filtering
#SBATCH --output=filtering_%A_%a.out
#SBATCH --error=filtering_%A_%a.err
#SBATCH -p shared
#SBATCH -n 6
#SBATCH -N 1
#SBATCH --time=1-00:00:00 
#SBATCH --mail-type=all
#SBATCH --mail-user=sferreira.mafalda@gmail.com

ml load bioinfo-tools bcftools

WD="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data"
INPUTWD=${WD}"/00_unfiltered_vcfs"
METADWD=${WD}"/metadata"

cd ${INPUTWD}

bcftools view --samples-file ${METADWD}/samples.txt -O z -o chr1.27indv.vcf.gz 1.vcf.gz 
bcftools view --samples-file ${METADWD}/samples.txt -O z -o chr7.27indv.vcf.gz 7.vcf.gz 

bcftools index chr1.27indv.vcf.gz
bcftools index chr7.27indv.vcf.gz

bcftools stats -s - chr1.27indv.vcf.gz > chr1.27indv.vcf.stats
bcftools stats -s - chr7.27indv.vcf.gz > chr7.27indv.vcf.stats

~~~

The files above are filtered like this:

##FILTER=<ID=SnpGap,Description="SNP within 5 bp of indel,other">
##bcftools_filterVersion=1.21-75-ga35ad84f+htslib-1.21-29-g243e97ec
##bcftools_filterCommand=filter -Ou -g 5:indel,other AllMouseAUTO.vcf.gz; Date=Thu Apr 10 00:28:53 2025
##bcftools_viewVersion=1.21-75-ga35ad84f+htslib-1.21-29-g243e97ec
##bcftools_viewCommand=view -Oz -M 2 -V indels,other,mnps; Date=Thu Apr 10 00:28:53 2025
##bcftools_filterCommand=filter -Oz -e INFO/DP>5031 AllMouseAUTO_SNP_ONLY.vcf.gz; Date=Thu Apr 10 18:08:45 2025
##bcftools_filterCommand=filter -Oz -e MQ<30 AllMouseAUTO_SNP_ONLY_depthFilter5031.vcf.gz; Date=Fri Apr 11 16:55:52 2025
##bcftools_filterCommand=filter -Oz -e QUAL<30 AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30.vcf.gz; Date=Sat Apr 12 13:28:39 2025
##bcftools_filterCommand=filter -Oz -e 'MEAN(FMT/SP) >3' AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30.vcf.gz; Date=Sat Apr 12 17:08:29 2025
##bcftools_filterCommand=filter -Ou -S . -e FMT/GQ<20 AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3.vcf.gz; Date=Sun Apr 13 15:26:08 2025
##bcftools_filterCommand=filter -Oz -S . -e FMT/DP<5 -; Date=Sun Apr 13 15:26:08 2025

So, this is a basic filtering. For our PCA and faststrcuture stuff we should further filter:

by maf
by 1 SNP every 1Kb? 
by missing data

~~~bash

WD="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data"


# FILTER 1: miss 0.2

# Here we are excluding sites where missing data is higher than 0.2. Which means we must have a genotype call for at least 80% of the individuals in every case. This is quite stringent.
mkdir filter_minmiss20

bcftools view -e "F_MISSING > 0.2" -O z -o ${WD}/filter_minmiss20/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.vcf.gz AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.vcf.gz

bcftools index ${WD}/filter_minmiss20/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.vcf.gz

bcftools stats -s - ${WD}/filter_minmiss20/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.vcf.gz > ${WD}/filter_minmiss20/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.stats

# BIALLELIC SNPs
# FILTER 2: miss 0.2 and maf 0.05 
# Supposedly we have 27 samples, 54 haplotypes. So the minor allele frequency of 3 variants is 0.05 (3/54) or (2/54) so 0.03. We should be good in this ballpark

mkdir filter_minmiss20_biallelic_snps_maf5

# -m 2 option makes sure each site has at least 2 alleles
# -M 2 options makes sure each site has at most 2 alleles
# --min-af 0.05:minor each variant site should have a minor allele frequency of 0.05 or more (we exclude low frequency variants)
# -v snps keeps only snp variantion and excludes indels

bcftools view -m 2 -M 2 -v snps --min-af 0.05:minor -O z -o ${WD}/filter_minmiss20_biallelic_snps_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.vcf.gz ${WD}/filter_minmiss20/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.vcf.gz 

bcftools index ${WD}/filter_minmiss20_biallelic_snps_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.vcf.gz

bcftools stats -s - ${WD}/filter_minmiss20_biallelic_snps_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.vcf.gz ${WD}/filter_minmiss20_biallelic_snps_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.stats
~~~


Wait, maybe let's do this in vcftools:



~~~bash
#!/bin/bash
#SBATCH -A naiss2025-22-172
#SBATCH --job-name=vcftools
#SBATCH --output=logs/vcftools_%A_%a.out
#SBATCH --error=logs/vcftools_%A_%a.err
#SBATCH -p shared
#SBATCH --array=1,7
#SBATCH -n 6
#SBATCH -N 1
#SBATCH --time=1-00:00:00 
#SBATCH --mail-type=all
#SBATCH --mail-user=sferreira.mafalda@gmail.com

# Load:
ml load bioinfo-tools vcftools bcftools

# Which chromosome:
ChrName=chr${SLURM_ARRAY_TASK_ID}

# Working directories:
WD="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data"
cd ${WD}

input_vcf="01_filtered_vcfs/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader."${ChrName}".27indv.vcf.gz"


mkdir 02_filtered_vcfs_for_analysis


# FILTER 1: miss 0.2

# Here we are excluding sites where missing data is higher than 0.2. Which means we must have a genotype call for at least 80% of the individuals in every case. This is quite stringent.

#Here, ‘––out’ specifies the name/location for the log file. If no name is specified, the log file will be named “out.log”. With‘––recode’ it is specified that a new VCF should be generated. To keep all meta information in the genotype data we should specify ‘––recode-INFO-all’. The output is directed to standard out (stdout) using ‘––stdout’ and can thus be piped to another program or directly written to a new file. With ‘bgzip’ the piped output is compressed.

# Invariant sites:
echo "Generating invariant files"
vcftools --gzvcf ${input_vcf} --max-missing 0.2 --max-non-ref-af 0 --recode --recode-INFO-all --stdout | bgzip -c > 02_filtered_vcfs_for_analysis/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.inv."${ChrName}".27indv.vcf.gz

# Variant sites:
echo "Generating variant files"
# FILTER 2: BIALLELIC SNPS miss 0.2, maf 0.05

vcftools --gzvcf ${input_vcf} --max-missing 0.2 --min-alleles 2 --max-alleles 2 --maf 0.05 --recode --recode-INFO-all --stdout | bgzip -c > 02_filtered_vcfs_for_analysis/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf5.var."${ChrName}".27indv.vcf.gz

# FILTER 3: BIALLELIC SNPS miss 0.2, maf 0.03
vcftools --gzvcf ${input_vcf} --max-missing 0.2 --min-alleles 2 --max-alleles 2 --maf 0.05 --recode --recode-INFO-all --stdout | bgzip -c > 02_filtered_vcfs_for_analysis/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf3.var."${ChrName}".27indv.vcf.gz

# JOIN VARIANT and INVARIANT:
echo "Indexing resulting vcf files"

# Index files:
bcftools index 02_filtered_vcfs_for_analysis/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.inv."${ChrName}".27indv.vcf.gz
bcftools index 02_filtered_vcfs_for_analysis/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf5.var."${ChrName}".27indv.vcf.gz
bcftools index 02_filtered_vcfs_for_analysis/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf3.var."${ChrName}".27indv.vcf.gz

# Concatenate files:
bcftools concat --allow-overlaps 02_filtered_vcfs_for_analysis/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.inv."${ChrName}".27indv.vcf.gz 02_filtered_vcfs_for_analysis/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf3.var."${ChrName}".27indv.vcf.gz -O z -o AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf3.AllSites."${ChrName}".27indv.vcf.gz

# Concatenate files:
bcftools concat --allow-overlaps 02_filtered_vcfs_for_analysis/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.inv."${ChrName}".27indv.vcf.gz 02_filtered_vcfs_for_analysis/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf5.var."${ChrName}".27indv.vcf.gz -O z -o AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf5.AllSites."${ChrName}".27indv.vcf.gz

~~~


The above worked very well, though the invariatn sites are not very realistic since the initial vcf file only had snps... We need to regenerate the Allsites. 

But I will start working with these files for Fst, PCA and others.

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

# Load:
ml load bioinfo-tools vcftools plink

# 
WD="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data"
cd ${WD}

INPUTWD=${WD}"/02_filtered_vcfs_for_analysis"
OUTPUTWD=${WD}"/03_PCA/01_inputs"
OUTPUTPCAWD=${WD}"/03_PCA/02_results"

# If the directories don't exist, create them:
mkdir -p $OUTPUTWD
mkdir -p $OUTPUTPCAWD

# To do population structure, it is generally better to run a window and select 1 SNP every ~10kb. I will try other windows anyway?
# Generate input files:
vcftools --thin 10000 --gzvcf ${INPUTWD}/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf3.var.chr1.27indv.vcf.gz --stdout --recode --recode-INFO-all | bgzip -c > $OUTPUTWD/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf3.var.10kbSNPs.chr1.27indv.vcf.gz

vcftools --thin 10000 --gzvcf ${INPUTWD}/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf5.var.chr1.27indv.vcf.gz --stdout --recode --recode-INFO-all | bgzip -c > $OUTPUTWD/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf5.var.10kbSNPs.chr1.27indv.vcf.gz

# PLINK PCA:
for vcf in ${OUTPUTWD}/*vcf.gz;
do
plink --vcf "$vcf" --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca \
--out ${OUTPUTPCAWD}/"$(basename "$vcf" .vcf.gz)";
done; 

~~~

Download eigenvec and eigenval files:


scp -i ~/.ssh/id-ed25519-pdc mafaldaf@dardel.pdc.kth.se:/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data/03_PCA/02_results/\* ./

REsults are the same with maf3 or maf5. let's go for maf5!

Let's run faststructure

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

# Load:
ml load bioinfo-tools fastStructure plink

# Set working environment variables:
WD="/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data"
cd ${WD}

# We can use the same inputs the .bed .bim. .fam files generated in the pca step:
INPUTWD=${WD}"/03_PCA/01_inputs"
OUTPUTPCAWD=${WD}"/04_FastStructure/02_results"

# If the directories don't exist, create them:
mkdir -p $OUTPUTPCAWD

# FastStructure accepts the classical plink input files, which we can generate with plink itself:

plink --vcf ${INPUTWD}/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf5.var.10kbSNPs.chr1.27indv.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --out ${INPUTWD}/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf5.var.10kbSNPs.chr1.27indv

structure.py -K 2 --input=${INPUTWD}/AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3_soft_filtered_20_5_reheader.miss20.biallelic.maf5.var.10kbSNPs.chr1.27indv --output=${OUTPUTPCAWD}/K2.chr1 --format=bed 
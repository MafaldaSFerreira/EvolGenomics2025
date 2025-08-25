# 2025-08-18

Over the next two days, we will conduct a population structure analysis and examine summary statistics for our populations 


# Filtering vcf files

We are working with a vcf file for chromosome 1, containing called genotypes for 27 individuals. It is good practice to apply filters to the called genotypes before continuing with formal analysis. We will different levels of confidence for the accuracy of the genotypes that were called depending on depth of coverage (how many reads map to a particular position), genotype quality assigned by the software we used to call genotypes, if variants are SNPs or INDELs, etc.

So let's perform some general filtering. For this, we will use the program bcftools. Open your terminal and log-in to the server. The vcf file we are going to filter is the one you have generated already in the previous classes. If you can find it, it is placed on the following path:



~~~bash
vcf=/cfs/klemming/scratch/m/mafaldaf/Teaching/20250900_EvolGenomics/mouse_data/00_unfiltered_vcfs/subset_vcf

chr1.27indv.vcf.gz
chr1.27indv.vcf.gz.csi

~~~

We will first mark variants that have

##bcftools_viewCommand=view -Oz -M 2 -V indels,other,mnps; Date=Thu Apr 10 00:28:53 2025
##bcftools_filterCommand=filter -Oz -e INFO/DP>5031 AllMouseAUTO_SNP_ONLY.vcf.gz; Date=Thu Apr 10 18:08:45 2025
##bcftools_filterCommand=filter -Oz -e MQ<30 AllMouseAUTO_SNP_ONLY_depthFilter5031.vcf.gz; Date=Fri Apr 11 16:55:52 2025
##bcftools_filterCommand=filter -Oz -e QUAL<30 AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30.vcf.gz; Date=Sat Apr 12 13:28:39 2025
##bcftools_filterCommand=filter -Oz -e 'MEAN(FMT/SP) >3' AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30.vcf.gz; Date=Sat Apr 12 17:08:29 2025
##bcftools_filterCommand=filter -Ou -S . -e FMT/GQ<20 AllMouseAUTO_SNP_ONLY_depthFilter5031_MQ30_QUAL30_SP_3.vcf.gz; Date=Sun Apr 13 15:26:08 2025
##bcftools_filterCommand=filter -Oz -S . -e FMT/DP<5 -; Date=Sun Apr 13 15:26:08 2025



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

#bcftools filter -i 'QUAL>=30 && MQ>=30 && INFO/DP<1900' -Oz -o chr1.27indv.QUAL30_MQ30_depthFilter1900.vcf.gz chr1.27indv.vcf.gz

# Sample specific filters:
bcftools filter -S . -i 'TYPE="snp" & FMT/GQ>=20 & FMT/DP>5' -Oz -o chr1.27indv.QUAL30_MQ30_depthFilter1900_GQ20_DP5.vcf.gz chr1.27indv.QUAL30_MQ30_depthFilter1900.vcf.gz

# Finally, select variants that have PASS the above filters, are not indels, don't have * ALT alleles and where the reference is not a missing genotype "N"
bcftools view -f PASS -e 'ALT="*" | TYPE~"indel" | ref="N"' -O z -o chr1.27indvs.QUALFilters.vcf.gz chr1.27indv.QUAL30_MQ30_depthFilter1900_GQ20_DP5.vcf.gz

bcftools index chr1.27indvs.QUALFilters.vcf.gz
bcftools stats -s - chr1.27indvs.QUALFilters.vcf.gz > chr1.27indvs.QUALFilters.stats
~~~


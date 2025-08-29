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

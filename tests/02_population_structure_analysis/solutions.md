## Genetic differentiation between populations and levels of diveristy 

Hello! I can see you are stuck in calculating pairwise Fst between populations using `vcftools`! Don't worry. Look at my solution below. If you are working in Duke, check the first solution.

We will  use `vcftools` to make the calculations we need. We will perform calculations in a pairwise fashion between populations. We will use several input files contained in `metadata_files/` which contain a list of individual names per population.


### For Duke 


~~~bash
cd /data/MsC2025/maffer/population_analysis
~~~

~~~bash
cp /data/MsC2025/shared_data/practical_4/metadata_files/GER.txt .
cp /data/MsC2025/shared_data/practical_4/metadata_files/FRA.txt .
cp /data/MsC2025/shared_data/practical_4/metadata_files/IND.txt .
cp /data/MsC2025/shared_data/practical_4/metadata_files/HEL.txt .
~~~

Between GER vs HEL:
~~~bash
vcftools --gzvcf chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.vcf.gz --weir-fst-pop GER.txt --weir-fst-pop HEL.txt --out GER_vs_HEL
~~~

The remaining comparisons:

~~~bash 
vcftools --gzvcf chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.vcf.gz --weir-fst-pop GER.txt --weir-fst-pop FRA.txt --out GER_vs_FRA

vcftools --gzvcf chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.vcf.gz --weir-fst-pop GER.txt --weir-fst-pop IND.txt --out GER_vs_IND

vcftools --gzvcf chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.vcf.gz --weir-fst-pop FRA.txt --weir-fst-pop HEL.txt --out FRA_vs_HEL

vcftools --gzvcf chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.vcf.gz --weir-fst-pop FRA.txt --weir-fst-pop IND.txt --out FRA_vs_IND

vcftools --gzvcf chr1.27indvs.QUALFilters.POPFilters.10kbSNPs.vcf.gz --weir-fst-pop IND.txt --weir-fst-pop HEL.txt --out IND_vs_HEL
~~~

After the commands have run, you will find a set of `.log` files which will contain the Fst calculations per population.

Use `grep` to extract them all comparisons at the same time.

~~~
grep "Weir and Cockerham mean Fst estimate:" *log | cut -f1,3 -d":" | sed 's/.log:/\t/' > meanFst.out
grep "Weir and Cockerham weighted Fst estimate:" *log | cut -f1,3 -d":" | sed 's/.log:/\t/' > wFst.out
~~~

Use `cat` to inspect the file.



### For Dardel
~~~bash
ls metadata_files/
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
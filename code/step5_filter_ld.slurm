#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mac9jc@virginia.edu
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6000
#SBATCH --account=tumi
#SBATCH --output=/home/mac9jc/cparvum_genomes/outfiles/step5_ld.out

module purge
module load gatk/4.0.0.0

## rivanna paths
analyze_genome_path="/scratch/mac9jc/genome_analysis_project/results"
cd $analyze_genome_path
reference=/home/mac9jc/cparvum_genomes/genomes_reference/Cp_ref.fa

# SNPs
gatk SelectVariants -V output_ld.vcf -select-type SNP -O snps_ld.vcf.gz

# FILTER
gatk VariantFiltration \
-V snps_ld.vcf.gz \
-filter "QD < 15.0" --filter-name "QD15" \
-filter "FS > 12.0" --filter-name "FS12" \
-filter "MQ < 58.0" --filter-name "MQ58" \
-filter "MQRankSum < -3.0" --filter-name "MQRankSum-3" \
-filter "ReadPosRankSum < -3.0" --filter-name "ReadPosRankSum-3" \
-O all_ld_filtered.vcf.gz

gatk SelectVariants \
-R $reference \
-V all_ld_filtered.vcf.gz \
--exclude-filtered true	\
-O all_ld_filtered_removed.vcf.gz


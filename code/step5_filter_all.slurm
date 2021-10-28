#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mac9jc@virginia.edu
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6000
#SBATCH --account=tumi
#SBATCH --output=/home/mac9jc/cparvum_genomes/outfiles/step5_all.out

module purge
module load gatk/4.0.0.0

## rivanna paths
analyze_genome_path="/scratch/mac9jc/genome_analysis_project/results"
cd $analyze_genome_path
reference=/home/mac9jc/cparvum_genomes/genomes_reference/Cp_ref.fa

# SNPs
gatk SelectVariants -V output.vcf -select-type SNP -O snps.vcf.gz

# INDELs
gatk SelectVariants -V output.vcf -select-type INDEL -O indels.vcf.gz

# FILTER
gatk VariantFiltration \
-V snps.vcf.gz \
-filter "QUAL < 25.0"  --filter-name "Qual" \
-filter "QD < 15.0" --filter-name "QD15" \
-filter "FS > 12.0" --filter-name "FS12" \
-filter "MQ < 58.0" --filter-name "MQ58" \
-filter "MQRankSum < -3.0" --filter-name "MQRankSum-3" \
-filter "ReadPosRankSum < -3.0" --filter-name "ReadPosRankSum-3" \
-O snps_filtered.vcf.gz

gatk SelectVariants \
-R $reference \
-V snps_filtered.vcf.gz \
--exclude-filtered true \
-O snps_filtered_removed.vcf.gz

gatk VariantFiltration \
-V indels.vcf.gz \
-filter "QUAL < 25.0"  --filter-name "Qual" \
-filter "QD < 15.0" --filter-name "QD15" \
-filter "FS > 12.0" --filter-name "FS12" \
-filter "MQ < 58.0" --filter-name "MQ58" \
-filter "MQRankSum < -3.0" --filter-name "MQRankSum-3" \
-filter "ReadPosRankSum < -3.0" --filter-name "ReadPosRankSum-3" \
-O indels_filtered.vcf.gz

gatk SelectVariants \
-R $reference \
-V indels_filtered.vcf.gz \
--exclude-filtered true \
-O indels_filtered_removed.vcf.gz


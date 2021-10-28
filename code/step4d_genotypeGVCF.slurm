#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mac9jc@virginia.edu
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=12000
#SBATCH --account=tumi
#SBATCH --output=/home/mac9jc/cparvum_genomes/outfiles/step4_p4.out

module purge
module load gatk/4.0.0.0

## rivanna paths
analyze_genome_path="/scratch/mac9jc/genome_analysis_project/results"
ref_genome_path="/home/mac9jc/cparvum_genomes/genomes_reference"
cd $analyze_genome_path

gatk --java-options "-Xmx62G" GenotypeGVCFs \
   --reference $ref_genome_path/Cp_ref.fa \
   --variant /scratch/mac9jc/genome_analysis_project/results/all.vcf \
   --output /scratch/mac9jc/genome_analysis_project/results/output.vcf \
   -new-qual true

gatk --java-options "-Xmx62G" GenotypeGVCFs \
   --reference $ref_genome_path/Cp_ref.fa \
   --variant /scratch/mac9jc/genome_analysis_project/results/all_ld.vcf \
   --output /scratch/mac9jc/genome_analysis_project/results/output_ld.vcf \
   -new-qual true


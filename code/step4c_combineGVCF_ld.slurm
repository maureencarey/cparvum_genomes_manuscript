#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mac9jc@virginia.edu
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6000
#SBATCH --account=tumi
#SBATCH --output=/home/mac9jc/cparvum_genomes/outfiles/step4_p3_ld.out

module purge
module load gatk/4.0.0.0

## rivanna paths
analyze_genome_path="/scratch/mac9jc/genome_analysis_project/results"
ref_genome_path="/home/mac9jc/cparvum_genomes/genomes_reference"
cd $analyze_genome_path

gatk --java-options "-Xmx4G" CombineGVCFs --reference $ref_genome_path/Cp_ref.fa \
   --variant /scratch/mac9jc/genome_analysis_project/results/icddrb29_RG_ld.g.vcf \
   --variant /scratch/mac9jc/genome_analysis_project/results/icddrb47_RG_ld.g.vcf \
   --variant /scratch/mac9jc/genome_analysis_project/results/icddrb63_RG_ld.g.vcf \
   --variant /scratch/mac9jc/genome_analysis_project/results/icddrb90_RG_ld.g.vcf \
   --variant /scratch/mac9jc/genome_analysis_project/results/icddrb93_RG_ld.g.vcf \
   --variant /scratch/mac9jc/genome_analysis_project/results/icddrb111_RG_ld.g.vcf \
   --variant /scratch/mac9jc/genome_analysis_project/results/SRR6147587_RG_ld.g.vcf \
   --variant /scratch/mac9jc/genome_analysis_project/results/SRR6147945_RG_ld.g.vcf \
   --variant /scratch/mac9jc/genome_analysis_project/results/SRR6147964_RG_ld.g.vcf \
   --variant /scratch/mac9jc/genome_analysis_project/results/SRR6148259_RG_ld.g.vcf \
   --variant /scratch/mac9jc/genome_analysis_project/results/SRR6813716_RG_ld.g.vcf \
   --variant /scratch/mac9jc/genome_analysis_project/results/SRR6813717_RG_ld.g.vcf \
   --variant /scratch/mac9jc/genome_analysis_project/results/SRR6813718_RG_ld.g.vcf \
   --variant /scratch/mac9jc/genome_analysis_project/results/SRR6813719_RG_ld.g.vcf \
   --variant /scratch/mac9jc/genome_analysis_project/results/SRR7898459_RG_ld.g.vcf \
   --variant /scratch/mac9jc/genome_analysis_project/results/SRR6117460_RG_ld.g.vcf \
   --variant /scratch/mac9jc/genome_analysis_project/results/SRR6147472_RG_ld.g.vcf \
   --variant /scratch/mac9jc/genome_analysis_project/results/SRR6147581_RG_ld.g.vcf \
   --output /scratch/mac9jc/genome_analysis_project/results/all_ld.vcf

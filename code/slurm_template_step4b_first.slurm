#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mac9jc@virginia.edu
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=4
#SBATCH --account=tumi

module purge
module load gatk/4.0.0.0
module load samtools

## rivanna paths
analyze_genome_path="/scratch/mac9jc/genome_analysis_project/results"
ref_genome_path="/home/mac9jc/cparvum_genomes/genomes_reference"
cd $analyze_genome_path


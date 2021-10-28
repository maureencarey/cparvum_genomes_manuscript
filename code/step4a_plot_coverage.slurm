#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mac9jc@virginia.edu
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64000
#SBATCH --account=tumi
#SBATCH --output=/home/mac9jc/cparvum_genomes/outfiles/plot_coverage.out

module purge
module load gcc/7.1.0 openmpi/3.1.4 intel/18.0  intelmpi/18.0 R/4.0.0
module load samtools

## rivanna paths
analyze_genome_path="/scratch/mac9jc/genome_analysis_project/results"
cd $analyze_genome_path

FILES=$analyze_genome_path/*_dedup.bam
for f in $FILES; do 
    output=${f%.*}"coverage.txt"
    samtools depth $f > $output
done

Rscript /home/mac9jc/cparvum_genomes/code/make_coverage_plot.R


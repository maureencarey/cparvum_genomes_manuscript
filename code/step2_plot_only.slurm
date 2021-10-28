#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mac9jc@virginia.edu
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=1
#SBATCH --account=tumi
#SBATCH --output=/home/mac9jc/cparvum_genomes/outfiles/step2_prealignmentQC_plot_only.out

module purge
module load gcc/7.1.0  openmpi/3.1.4 intel/18.0  intelmpi/18.0
module load R/4.0.0
Rscript ~/cparvum_genomes/analyses/get_input_reads.R




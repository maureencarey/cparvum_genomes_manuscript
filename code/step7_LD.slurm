#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mac9jc@virginia.edu
#SBATCH --partition=largemem
#SBATCH --mem-per-cpu=62000
#SBATCH --cpus-per-task=1
#SBATCH --account=tumi
#SBATCH --output=/home/mac9jc/cparvum_genomes/outfiles/step7_LD.out

module purge
module load gatk/4.0.0.0
module load plink/1.90b6.16
module load bcftools
module load vcftools/0.1.16

bash /home/mac9jc/cparvum_genomes/analyses/LD.sh

module purge
module load gcc/7.1.0  openmpi/3.1.4
module load R/4.0.0

Rscript /home/mac9jc/cparvum_genomes/analyses/plotLD.R 

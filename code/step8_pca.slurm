#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --partition=largemem
#SBATCH --mem-per-cpu=62000
#SBATCH --cpus-per-task=1
#SBATCH --account=tumi
#SBATCH --output=temp_plot_pca.out

module purge
module load plink/1.90b6.16

plink --vcf /scratch/mac9jc/genome_analysis_project/results/all_ld_filtered_removed.vcf.gz --out pca \
--allow-no-sex --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --chr-set 8 no-y no-xy no-mt -mind 0.2 --pca

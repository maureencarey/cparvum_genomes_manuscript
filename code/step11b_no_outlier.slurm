#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --partition=largemem
#SBATCH --mem-per-cpu=62000
#SBATCH --cpus-per-task=1
#SBATCH --account=tumi
#SBATCH --output=/home/mac9jc/cparvum_genomes/outfiles/step11b.out


module purge
module load plink/1.90b6.16
cd /home/mac9jc/cparvum_genomes/pca_results

# manually remove rows with gp60 cgd6_1080
plink --vcf /scratch/mac9jc/genome_analysis_project/results/no_47_snps_updated_chrom_annotated_missense_secreted_copy_no_gp60.vcf --out no_47_pca_missense_secreted_no_gp60 \
--allow-no-sex --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --chr-set 8 no-y no-xy no-mt -mind 0.2 --pca

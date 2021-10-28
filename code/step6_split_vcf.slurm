#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mac9jc@virginia.edu
#SBATCH --partition=dev
#SBATCH --cpus-per-task=1
#SBATCH --account=tumi
#SBATCH --output=/home/mac9jc/cparvum_genomes/outfiles/step6_split.out

module purge
module load vcftools/0.1.16
module load bcftools/1.9

sra_genome_path="/scratch/mac9jc/genome_analysis_project/results"
cd $sra_genome_path

# make individual sample files
file='/scratch/mac9jc/genome_analysis_project/results/snps_filtered_removed.vcf.gz'
for sample in `bcftools query -l $file`; do
    echo $sample
    bcftools view -c1 -Oz -s $sample -o ${file/.vcf*/.$sample.vcf.gz} $file
done

# make individual sample group files
# Bangladesh does not include mixed sample
cd /scratch/mac9jc/genome_analysis_project/results/
FILES1=/scratch/mac9jc/genome_analysis_project/results/snps_filtered_removed.*.vcf.gz
for f in $FILES1; do
    bcftools index $f -f
done
# doesn't work with polyploid
#bcftools merge --missing-to-ref 

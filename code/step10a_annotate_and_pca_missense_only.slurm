#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --partition=largemem
#SBATCH --mem-per-cpu=62000
#SBATCH --cpus-per-task=1
#SBATCH --account=tumi
#SBATCH --output=/home/mac9jc/cparvum_genomes/outfiles/step10_missense.out

module purge
module load bcftools/1.9

echo 'first rename chromosomes'

# rename chromosomes to be compatible with genome reference (below)
cd /scratch/mac9jc/genome_analysis_project/results/
bcftools annotate --rename-chrs /home/mac9jc/cparvum_genomes/chromosome_dict2.txt all_ld_filtered_removed.vcf.gz -o snps_updated_chrom.vcf.gz

echo 'starting snp annotation'

# annotate all samples
cp snps_updated_chrom.vcf.gz snps_updated_chrom.vcf
java -jar /project/tumi/carey/snpEff/snpEff/snpEff.jar Cryptosporidium_parvum_iowa_ii_gca_000165345 snps_updated_chrom.vcf > snps_updated_chrom_annotated.vcf

echo 'starting snp extraction'

##### get missense only ######

# extract missense mutations only
java -jar /project/tumi/carey/snpEff/snpEff/SnpSift.jar filter "ANN[*].EFFECT has 'missense_variant'" snps_updated_chrom_annotated.vcf > snps_updated_chrom_annotated_missense.vcf

# get number of SNPs
cat snps_updated_chrom_annotated.vcf | grep -v ^# | wc -l
cat snps_updated_chrom_annotated_missense.vcf | grep -v ^# | wc -l

##### get missense in secreted proteins only ######

# extract missense mutations in list of genes
java -jar /project/tumi/carey/snpEff/snpEff/SnpSift.jar filter -s ~/cparvum_genomes/secreted_cparvum_genes_ids.csv "((ANN[*].EFFECT has 'missense_variant') & (ANN[*].GENEID in SET[0]))" snps_updated_chrom_annotated.vcf > snps_updated_chrom_annotated_missense_secreted.vcf

#####

module purge
module load plink/1.90b6.16

echo 'starting pca'
cd /home/mac9jc/cparvum_genomes/pca_results

plink --vcf /scratch/mac9jc/genome_analysis_project/results/snps_updated_chrom_annotated_missense.vcf --out pca_missense \
--allow-no-sex --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --chr-set 8 no-y no-xy no-mt -mind 0.2 --pca

echo 'starting second pca'

plink --vcf /scratch/mac9jc/genome_analysis_project/results/snps_updated_chrom_annotated_missense_secreted.vcf --out pca_missense_secreted \
--allow-no-sex --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --chr-set 8 no-y no-xy no-mt -mind 0.2 --pca

echo 'starting last pca'
plink --vcf /scratch/mac9jc/genome_analysis_project/results/all_ld_filtered_removed.vcf.gz --out pca \
--allow-no-sex --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --chr-set 8 no-y no-xy no-mt -mind 0.2 --pca

echo 'starting final (for real) pca'
cp /scratch/mac9jc/genome_analysis_project/results/snps_updated_chrom_annotated_missense_secreted.vcf /scratch/mac9jc/genome_analysis_project/results/snps_updated_chrom_annotated_missense_secreted_copy_no_gp60.vcf
# manually remove rows with gp60 cgd6_1080 

#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --partition=largemem
#SBATCH --mem-per-cpu=62000
#SBATCH --cpus-per-task=1
#SBATCH --account=tumi
#SBATCH --output=/home/mac9jc/cparvum_genomes/outfiles/step11a.out


# this runs step 4c-10 on all sampels except icddrb47 since it is an outlier in terms of number of unique snps (high), low read count, and via the PCA
module purge


# step 4c

module load gatk/4.0.0.0

## rivanna paths
analyze_genome_path="/scratch/mac9jc/genome_analysis_project/results"
ref_genome_path="/home/mac9jc/cparvum_genomes/genomes_reference"
cd $analyze_genome_path

gatk --java-options "-Xmx4G" CombineGVCFs --reference $ref_genome_path/Cp_ref.fa \
   --variant /scratch/mac9jc/genome_analysis_project/results/icddrb29_RG_ld.g.vcf \
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
   --output /scratch/mac9jc/genome_analysis_project/results/all_ld_no_47.vcf

# step 4d

analyze_genome_path="/scratch/mac9jc/genome_analysis_project/results"
ref_genome_path="/home/mac9jc/cparvum_genomes/genomes_reference"
cd $analyze_genome_path

gatk --java-options "-Xmx62G" GenotypeGVCFs \
   --reference $ref_genome_path/Cp_ref.fa \
   --variant /scratch/mac9jc/genome_analysis_project/results/all_ld_no_47.vcf \
   --output /scratch/mac9jc/genome_analysis_project/results/output_ld_no_47.vcf \
   -new-qual true

# step 5


analyze_genome_path="/scratch/mac9jc/genome_analysis_project/results"
cd $analyze_genome_path
reference=/home/mac9jc/cparvum_genomes/genomes_reference/Cp_ref.fa

# SNPs
gatk SelectVariants -V output_ld_no_47.vcf -select-type SNP -O no_47_snps_ld.vcf.gz

# FILTER
gatk VariantFiltration \
-V no_47_snps_ld.vcf.gz \
-filter "QD < 15.0" --filter-name "QD15" \
-filter "FS > 12.0" --filter-name "FS12" \
-filter "MQ < 58.0" --filter-name "MQ58" \
-filter "MQRankSum < -3.0" --filter-name "MQRankSum-3" \
-filter "ReadPosRankSum < -3.0" --filter-name "ReadPosRankSum-3" \
-O no_47_all_ld_filtered.vcf.gz

gatk SelectVariants \
-R $reference \
-V no_47_all_ld_filtered.vcf.gz \
--exclude-filtered true \
-O no_47_all_ld_filtered_removed.vcf.gz


# step 11

module load bcftools/1.9

echo 'first rename chromosomes'

# rename chromosomes to be compatible with genome reference (below)
cd /scratch/mac9jc/genome_analysis_project/results/
bcftools annotate --rename-chrs /home/mac9jc/cparvum_genomes/chromosome_dict2.txt no_47_all_ld_filtered_removed.vcf.gz -o no_47_snps_updated_chrom.vcf.gz

echo 'starting snp annotation'

# annotate all samples
cp no_47_snps_updated_chrom.vcf.gz no_47_snps_updated_chrom.vcf
java -jar /project/tumi/carey/snpEff/snpEff/snpEff.jar Cryptosporidium_parvum_iowa_ii_gca_000165345 no_47_snps_updated_chrom.vcf > no_47_snps_updated_chrom_annotated.vcf

echo 'starting snp extraction'

##### get missense only ######

# extract missense mutations only
java -jar /project/tumi/carey/snpEff/snpEff/SnpSift.jar filter "ANN[*].EFFECT has 'missense_variant'" no_47_snps_updated_chrom_annotated.vcf > no_47_snps_updated_chrom_annotated_missense.vcf

##### get missense in secreted proteins only ######

# extract missense mutations in list of genes
java -jar /project/tumi/carey/snpEff/snpEff/SnpSift.jar filter -s ~/cparvum_genomes/secreted_cparvum_genes_ids.csv "((ANN[*].EFFECT has 'missense_variant') & (ANN[*].GENEID in SET[0]))" no_47_snps_updated_chrom_annotated.vcf > no_47_snps_updated_chrom_annotated_missense_secreted.vcf

#####

module purge
module load plink/1.90b6.16

echo 'starting pca'
cd /home/mac9jc/cparvum_genomes/pca_results

plink --vcf /scratch/mac9jc/genome_analysis_project/results/no_47_snps_updated_chrom_annotated.vcf --out no_47_pca \
--allow-no-sex --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --chr-set 8 no-y no-xy no-mt -mind 0.2 --pca

cp /scratch/mac9jc/genome_analysis_project/results/no_47_snps_updated_chrom_annotated_missense_secreted.vcf /scratch/mac9jc/genome_analysis_project/results/no_47_snps_updated_chrom_annotated_missense_secreted_copy_no_gp60.vcf
# manually remove rows with gp60 cgd6_1080

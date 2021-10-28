## rivanna paths
analyze_genome_path="/scratch/mac9jc/genome_analysis_project/results"
cd $analyze_genome_path

vcftools --gzvcf all_ld_filtered_removed.vcf.gz --recode --recode-INFO-all --out temp_snps
gatk IndexFeatureFile -F temp_snps.recode.vcf
#plink --vcf temp_snps.recode.vcf --maf 0.05 --recode --out genotypes --allow-no-sex --double-id --allow-extra-chr
plink --vcf temp_snps.recode.vcf --recode --out genotypes --allow-no-sex --double-id --allow-extra-chr
# make ped file
plink --file genotypes --make-bed --out genotypes_ready --allow-no-sex --allow-extra-chr
# Fix bad SNP ids
Rscript /home/mac9jc/cparvum_genomes/code/makeSNPnamesUnique.R /scratch/mac9jc/genome_analysis_project/results/genotypes_ready.bim
# get allele frequency
plink --bfile genotypes_ready --freq --out freq_stat --allow-no-sex --chr-set 8 no-y no-xy no-mt --allow-extra-chr
# do LD
plink --bfile genotypes_ready --r2 --ld-window-r2 0.05 --ld-window 10 --ld-window-kb 30 --out LD_output --allow-no-sex --double-id --chr-set 8 no-y no-xy no-mt --allow-extra-chr

#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mac9jc@virginia.edu
#SBATCH --partition=dev
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6000
#SBATCH --account=tumi
#SBATCH --output=/home/mac9jc/cparvum_genomes/outfiles/step9_snpEff_all.out

module purge
module load bcftools/1.9

echo 'rename chrom'

# rename chromosomes to be compatible with genome reference (below)
cd /scratch/mac9jc/genome_analysis_project/results/
FILES=/scratch/mac9jc/genome_analysis_project/results/snps_filtered_removed.*.vcf.gz
for f in $FILES; do
   file_without_ext="${f%.*}"
   file_without_any_ext="${file_without_ext%.*}"
   renamed=$file_without_any_ext"_renamed.vcf.gz"
   bcftools annotate --rename-chrs /home/mac9jc/cparvum_genomes/chromosome_dict2.txt $f -o $renamed
done

echo 'get reference'

# make genome reference
cd /project/tumi/carey/snpEff/snpEff
java -jar snpEff.jar download -v Cryptosporidium_parvum_iowa_ii_gca_000165345

# annotate all samples' files
# http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf

echo 'starting to annotate'

cd /scratch/mac9jc/genome_analysis_project/results/
for f in $(ls snps_filtered_removed.*_renamed.vcf.gz); do
   # get rid of gz because not gzipped properly
   new_name="${f%.*}"
   mv $f $new_name
   output=$new_name".ann.vcf"
   # annotate with snpEff
   java -jar /project/tumi/carey/snpEff/snpEff/snpEff.jar Cryptosporidium_parvum_iowa_ii_gca_000165345 $new_name > $output
done

# extract missense mutations only
for f in $(ls snps_filtered_removed.*_renamed.vcf.ann.vcf); do
   output2=$f"_2.vcf"
   java -jar /project/tumi/carey/snpEff/snpEff/SnpSift.jar filter "ANN[*].EFFECT has 'missense_variant'" $f > $output2
   cp $output2 /scratch/mac9jc/genome_analysis_project/results/copy_to_mac/all
done

# get number of SNPs
#cat snps_filtered_removed.icddrb111.vcf.ann.vcf | grep -v ^# | wc -l
#cat snps_filtered_removed.icddrb111.vcf.ann.vcf_2.vcf | grep -v ^# | wc -l

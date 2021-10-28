#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mac9jc@virginia.edu
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6000
#SBATCH --account=tumi
#SBATCH --output=/home/mac9jc/cparvum_genomes/outfiles/step1_get_genomes.out

module purge
module load sratoolkit/2.9.1
module load gcc/7.1.0
module load fastqc/0.11.5

sra_genome_path="/scratch/mac9jc/genome_analysis_project"
cd $sra_genome_path

# get reference genomes and annotations
wget https://cryptodb.org/common/downloads/Current_Release/CparvumIowaII/fasta/data/CryptoDB-53_CparvumIowaII_Genome.fasta

module load gcc/9.2.0 bwa/0.7.17
module load samtools/1.9
module load picard/2.20.6

ref_genome_path="/home/mac9jc/cparvum_genomes/genomes_reference"
cd $ref_genome_path
mv $sra_genome_path/CryptoDB-53_CparvumIowaII_Genome.fasta Cp_ref.fa

# Make reference genome indexes for alignment steps
bwa index Cp_ref.fa
#java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R= Cp_ref.fa
samtools faidx Cp_ref.fa

# get publically available raw reads
cd /scratch/mac9jc/genome_analysis_project/
FILES="SRR6117460 SRR6147472 SRR6147581 SRR6147587 SRR6147945 SRR6147964 SRR6148259 SRR6813717 SRR6813716 SRR6813718 SRR6813719 SRR7898459"
#for f in $FILES; do fastq-dump --split-files --gzip $f; done
echo 'got files'

# get our genomes
#cp -R /home/mac9jc/cparvum_genomes/genomes_icddrb/. ./
#echo 'moved files'

fastqc_output_path="/scratch/mac9jc/genome_analysis_project/Cparvum_fastQC_output"
# evaluate quality of reads
#fastqc $sra_genome_path/*.fastq.gz 
echo 'ran fastqc'

#!/bin/bash

module purge
#module load gcc/7.1.0
#module load bwa/0.7.17
#module load samtools/1.9

ref_genome_path="/home/mac9jc/cparvum_genomes/genomes_reference"

# Make reference genome indexes for alignment steps
cd $ref_genome_path
#bwa index Cp_ref.fa
##java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R= Cp_ref.fa
#samtools faidx Cp_ref.fa
cd /home/mac9jc/cparvum_genomes/code

files="icddrb29 icddrb47 icddrb63 icddrb90 icddrb93 icddrb111 SRR6117460 SRR6147472 SRR6147581 SRR6147587 SRR6147945 SRR6147964 SRR6148259 SRR6813716 SRR6813717 SRR6813718 SRR6813719 SRR7898459"
#FILES=$sra_genome_path/*_1_dedup.fq

for filename in $files; do
    echo "$filename"

    #construct name for sbatch file being generated
    foo=${filename}"_step3.slurm"
    cp step3_template_p1.slurm $foo
    echo "#SBATCH --output=/scratch/mac9jc/genome_analysis_project/outfiles/"${filename}"_step3.out" >> $foo
    cat step3_template_p2.slurm >> $foo
#    echo "FILES="${filename}"_1_dedup.fq" >> $foo
    echo ${filename}"_1_dedup.fq" >> $foo
    cat step3_template_p3.slurm >> $foo
    sbatch $foo
done

#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mac9jc@virginia.edu
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=4
#SBATCH --account=tumi

module purge
module load picard/2.20.6

ref_genome_path="/home/mac9jc/cparvum_genomes/genomes_reference"
cd $ref_genome_path
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R= Cp_ref.fa

cd /home/mac9jc/cparvum_genomes/outfiles
for filename in /scratch/mac9jc/genome_analysis_project/results/*_RG.bam; do
    echo "$filename"
    file_without_pre="${filename##*/}"
    file_without_ext="${file_without_pre%.*}"
    species_string="${file_without_ext:0:${#file_without_ext}-3}"
    echo "$species_string"

    #construct name for sbatch file being generated
    foo="/home/mac9jc/cparvum_genomes/code/"${species_string}"_step4b.slurm"
    cp /home/mac9jc/cparvum_genomes/code/slurm_template_step4b_first.slurm $foo
    echo "f=$filename" >> $foo
    cat /home/mac9jc/cparvum_genomes/code/slurm_template_step4b_second.slurm >> $foo
    sbatch $foo
done

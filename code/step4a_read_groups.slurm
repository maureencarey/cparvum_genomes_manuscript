#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mac9jc@virginia.edu
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6000
#SBATCH --account=tumi
#SBATCH --output=/home/mac9jc/cparvum_genomes/outfiles/step4a_read_group.out

module purge
module load picard/2.20.6

## rivanna paths
analyze_genome_path="/scratch/mac9jc/genome_analysis_project/results"
ref_genome_path="/home/mac9jc/cparvum_genomes/genomes_reference"
cd $analyze_genome_path

FILES=$analyze_genome_path/*_dedup.bam
for f in $FILES; do
    file_without_ext="${f%.*}"
    f_without_path=$(basename $f)
    abbrev=${f_without_path%"_dedup.bam"}
    if [[ "$abbrev" =~ ^i.* ]]; then
        lib=$'lib1'
        if [[ "$abbrev" =~ ^icddrb9.* ]]; then
            use_for_RGPU=$'HX6'
        elif [[ "$abbrev" =~ ^icddrb1.* ]]; then
            use_for_RGPU=$'HX6'
        else
            use_for_RGPU=$'HX7'
        fi
    else
        lib='lib2'
        use_for_RGPU=$'unknown'
    fi
    java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=$f O=$abbrev"_RG.bam"\
        RGLB=$lib \
        RGPL=illumina \
        RGPU=$use_for_RGPU \
        RGSM=$abbrev
done

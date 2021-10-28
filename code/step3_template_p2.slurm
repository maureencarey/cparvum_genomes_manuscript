
module purge
module load gcc/9.2.0 
module load bwa/0.7.17
module load samtools/1.9
module load picard/2.20.6
module load bbmap/38.57

sra_genome_path="/scratch/mac9jc/genome_analysis_project/results"
ref_genome_path="/home/mac9jc/cparvum_genomes/genomes_reference"

# align
cd $sra_genome_path
FILES="

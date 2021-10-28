#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mac9jc@virginia.edu
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=1
#SBATCH --account=tumi
#SBATCH --output=/home/mac9jc/cparvum_genomes/outfiles/step2_prealignmentQC.out

module purge
module load gcc/9.2.0 bbmap/38.57
module load fastqc/0.11.5

sra_genome_path="/scratch/mac9jc/genome_analysis_project"
analyze_genome_path="/scratch/mac9jc/genome_analysis_project/results"
fastqc_output_path="/scratch/mac9jc/genome_analysis_project/Cparvum_fastQC_output"
cd $sra_genome_path

FILES1=$sra_genome_path/*.fastq.gz
FILES2=$sra_genome_path/*.fastq
for f in $FILES1 $FILES2; do
   echo '_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _'
   file_without_ext="${f%.*}"
   f_without_path=$(basename $f)
   # evaluate quality of reads
   fastqc $f

   if [[ $f == *_1.fast* ]]
   then
      echo "   trimming to 150 bp for: $f (well really trimming to multiple of 5 bp length"
      f_string=$(echo $(basename $file_without_ext)| cut -d'.' -f 1)
      f_without_num=${f_string%"_1"}
      f_rev=$f_without_num"_2.fastq.gz"
      output_file1=$analyze_genome_path/$f_without_num"_1_Ftm_trimmed.fq"
      output_file2=$analyze_genome_path/$f_without_num"_2_Ftm_trimmed.fq"
      bbduk.sh in1=$f in2=$f_rev out1=$output_file1 out2=$output_file2 ftm=5

      ## Adapter trimming # always recommended prior to merging
      echo "   adapter trimming for: $f"
      output_file3=$analyze_genome_path/$f_without_num"_1_Adapter_trimmed.fq"
      output_file4=$analyze_genome_path/$f_without_num"_2_Adapter_trimmed.fq"
      bbduk.sh in1=$output_file1 in2=$output_file2 out1=$output_file3 out2=$output_file4 ref=~/cparvum_genomes/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo

      ## Filter for spike ins LOOK UP SEQUENCING TECHNOLOGY TO SEE IF SOMETHING ELSE NEEDED, check names
      echo "   filter for spike ins for: $f"
      output_file5=$analyze_genome_path/$f_without_num"_1_Spikein_filtered.fq"
      output_file6=$analyze_genome_path/$f_without_num"_2_Spikein_filtered.fq"
      bbduk.sh in1=$output_file3 in2=$output_file4 out1=$output_file5 out2=$output_file6 ref=~/cparvum_genomes/phix_adapters.fa k=31 hdist=1 stats=$file_without_ext"_stats.txt"

      ## Quality filtering based on quality scores of end bases:  # sometimes not recommended prior to merging
      echo "   quality filtering based on end bases for: $f"
      output_file7=$analyze_genome_path/$f_without_num"_1_Q1.fq"
      output_file8=$analyze_genome_path/$f_without_num"_2_Q1.fq"
      bbduk.sh in1=$output_file5 in2=$output_file6 out1=$output_file7 out2=$output_file8 qtrim=rl trimq=10

      ## Quality filtering based on full read quality scores: # sometimes not recommended prior to merging
      echo "   quality filtering based on full read quality scores for: $f"
      output_file9=$analyze_genome_path/$f_without_num"_1_Q2.fq"
      output_file10=$analyze_genome_path/$f_without_num"_2_Q2.fq"
      bbduk.sh in1=$output_file7 in2=$output_file8 out1=$output_file9 out2=$output_file10 maq=10

      ## Additional quality filtering for SRA genomes - trim first 15 bp
      echo "   quality filtering based on end bases for: $f"
      output_file12=$analyze_genome_path/$f_without_num"_1_dedup.fq"
      output_file13=$analyze_genome_path/$f_without_num"_2_dedup.fq"
      if [ ${f_without_num:0:3} = "SRR" ]; then
          bbduk.sh in1=$output_file9 in2=$output_file10 out1=$output_file12 out2=$output_file13 ftl=15
      else
          mv $output_file9 $output_file12
          mv $output_file10 $output_file13
      fi

      # evaluate quality of reads
      fastqc $output_file12
      fastqc $output_file13

   else
      echo 'skip, rev reads'
   fi
done


mv *_fastqc.* $fastqc_output_path


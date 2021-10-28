"
for f in $FILES; do
    echo "   Aligning  $f file..."
    file_without_ext="${f%.*}"
    file_without_num=${file_without_ext%"_1_dedup"}
    f_rev=$file_without_num"_2_dedup.fq"
    output_file2=$file_without_num".sam"
    output_file_a=$file_without_num"_repaired1.fq"
    output_file_b=$file_without_num"_repaired2.fq"
    output_file_c=$file_without_num"_singletons.fq"
    output_file4=$file_without_num"_temp.bam"
    output_file5a=$file_without_num"_alignment_metrics1.txt"
    output_file5b=$file_without_num"_alignment_metrics2.txt"
    output_file6=$file_without_num"_dedup.bam"
    dups=$file_without_num"_dedup_metrics.txt"

    # repair
    repair.sh in1=$f in2=$f_rev out1=$output_file_a out2=$output_file_b outs=$output_file_c repair

    # align
    bwa mem -M $ref_genome_path/Cp_ref.fa $output_file_a $output_file_b > $output_file2
#    bwa mem -M $ref_genome_path/Cp_ref.fa $f $f_rev > $output_file2

    # convert and sort SAM (-> BAM)
    java -jar $EBROOTPICARD/picard.jar SortSam INPUT=$output_file2 OUTPUT=$output_file4 SORT_ORDER=coordinate

    # benchmark
    samtools flagstat $output_file4 > $output_file5a

    # mark duplicates
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
       I=$output_file4 \
       O=$output_file6 \
       M=$dups \
       REMOVE_DUPLICATES=true \
       REMOVE_SEQUENCING_DUPLICATES=true

    # benchmark
    samtools flagstat $output_file6 > $output_file5b

done


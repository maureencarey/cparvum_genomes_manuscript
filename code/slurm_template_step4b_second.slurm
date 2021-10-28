echo $f
file_without_ext="${f%.*}"
output_file0=$file_without_ext"_reord.bam"
# index
samtools sort -o $output_file0 $f
echo $output_file0
samtools index $output_file0

# call
gatk --java-options "-Xmx4G" HaplotypeCaller  \
   --reference $ref_genome_path/Cp_ref.fa \
   --input $output_file0 \
   --sample-ploidy 4 \
   -ERC GVCF \
   --output $file_without_ext'.g.vcf' \
   --max-reads-per-alignment-start 50 # default

# call for LD
gatk --java-options "-Xmx4G" HaplotypeCaller  \
   --reference $ref_genome_path/Cp_ref.fa \
   --input $output_file0 \
   --sample-ploidy 2 \
   -ERC GVCF \
   --output $file_without_ext'_ld.g.vcf' \
   --max-reads-per-alignment-start 50 # default



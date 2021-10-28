# cparvum_genomes

icddrb29 is a coinfection

waiting to push fasta files until certain that our consent covers this

# Contact

Feedback and questions to Maureen Carey - mac9jc [at] virginia [dot] edu

# code instructions

    ### run in sequence

    sbatch /home/mac9jc/cparvum_genomes/code/step1_get_genomes.slurm
    sbatch /home/mac9jc/cparvum_genomes/code/step2_oneloop.slurm
    sbatch /home/mac9jc/cparvum_genomes/code/step2_plot_only.slurm
    # rm /scratch/mac9jc/genome_analysis_project/results/*repaired*.fq # if a re run
    bash /home/mac9jc/cparvum_genomes/code/step3_align.sh
    sbatch /home/mac9jc/cparvum_genomes/code/step4a_read_groups.slurm
    sbatch /home/mac9jc/cparvum_genomes/code/step4a_plot_coverage.slurm
    sbatch /home/mac9jc/cparvum_genomes/code/step4b_generate_slurm_files.slurm

    ### these two can be run concurrently
    sbatch /home/mac9jc/cparvum_genomes/code/step4c_combineGVCF.slurm
    sbatch /home/mac9jc/cparvum_genomes/code/step4c_combineGVCF_ld.slurm

    ### back to in sequence
    sbatch /home/mac9jc/cparvum_genomes/code/step4d_genotypeGVCF.slurm

    ### these two can be concurrently
    sbatch /home/mac9jc/cparvum_genomes/code/step5_filter_all.slurm
    sbatch /home/mac9jc/cparvum_genomes/code/step5_filter_ld.slurm

    ### back to in sequence
    sbatch /home/mac9jc/cparvum_genomes/code/step6_split_vcf.slurm
    sbatch /home/mac9jc/cparvum_genomes/code/step7_LD.slurm
    sbatch /home/mac9jc/cparvum_genomes/code/step8_pca.slurm
    sbatch /home/mac9jc/cparvum_genomes/code/step9_snpEff.slurm 
    sbatch /home/mac9jc/cparvum_genomes/code/step10_annotate_and_pca_missense_only.slurm 

    ### run on mac stuff
    # copy  to mac: /scratch/mac9jc/genome_analysis_project/results/copy_to_mac/all/*
    # upsetR_figures.R
    # pca_variants.R

    #indel length
    module load vcftools/0.1.16
    vcftools --gzvcf /scratch/mac9jc/genome_analysis_project/results/indels_filtered_removed.vcf.gz --out sample --hist-indel-len
    tail sample.indel.hist

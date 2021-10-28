library("ShortRead")
library("dada2")
library("ggplot2")
library("dplyr")
sink('/home/mac9jc/cparvum_genomes/R_outfile_get_reads.txt')

# get raw sample reads 
path1 = "/scratch/mac9jc/genome_analysis_project"
fns1 <- list.files(path1)
fns1 = grep(".fastq.gz",fns1,value=TRUE)

### Load forward and reverse reads
fastqs1 <- fns1[grepl(".fastq.gz$", fns1)]
fastqs1 <- sort(fastqs1) # Sort ensures forward/reverse reads are in same order
fnFs1 <- fastqs1[grepl("_1", fastqs1)] # Just the forward read files
fnRs1 <- fastqs1[grepl("_2", fastqs1)] # Just the reverse read files
# Get sample names from the first part of the forward read filenames
sample.names1 <- sapply(strsplit(fnFs1, "_"), `[`, 1)
# Fully specify the path for the fnFs and fnRs
fnFs1 <- file.path(path1, fnFs1)
fnRs1 <- file.path(path1, fnRs1)

print(sample.names1)
df = data.frame(sample = sample.names1, forward = NA, reverse = NA)

for (i in 1:length(sample.names1)) {
  f_seq = getSequences(fnFs1[i]); r_seq = getSequences(fnRs1[i])
  df[i,2] = length(f_seq); df[i,3] = length(r_seq)} 
write.csv(df, "/home/mac9jc/cparvum_genomes/results/starting_seqs.csv")

df_mean = df %>% summarize(mean = mean(forward), median = median(forward))
p = ggplot(df) + geom_histogram(aes(x = forward)) + ylab('samples') + xlab('reads') +
  geom_vline(data = df_mean, aes(xintercept = mean), color = "red", linetype = "dashed") +
  theme(axis.text.x = element_text(size = 6)) + scale_x_log10()
ggsave(paste0("/home/mac9jc/cparvum_genomes/plots/read_counts_"
              , format(Sys.time(), "%Y_%m_%d")
              , ".png"), p, width = 8, height = 2, dpi = 600)


# get processed sample reads 
path1 = "/scratch/mac9jc/genome_analysis_project/results"
fns1 <- list.files(path1)
fns1 = grep("_dedup.fq",fns1,value=TRUE)

### Load forward and reverse reads
fastqs1 <- fns1[grepl("_dedup.fq$", fns1)]
fastqs1 <- sort(fastqs1) # Sort ensures forward/reverse reads are in same order
fnFs1 <- fastqs1[grepl("_1_dedup", fastqs1)] # Just the forward read files
fnRs1 <- fastqs1[grepl("_2_dedup", fastqs1)] # Just the reverse read files
# Get sample names from the first part of the forward read filenames
sample.names1 <- sapply(strsplit(fnFs1, "_"), `[`, 1)
# Fully specify the path for the fnFs and fnRs
fnFs1 <- file.path(path1, fnFs1)
fnRs1 <- file.path(path1, fnRs1)

print(sample.names1)
df = data.frame(sample = sample.names1, forward = NA, reverse = NA)

for (i in 1:length(sample.names1)) {
  f_seq = getSequences(fnFs1[i]); r_seq = getSequences(fnRs1[i])
  df[i,2] = length(f_seq); df[i,3] = length(r_seq)} 
write.csv(df, "/home/mac9jc/cparvum_genomes/results/processed_seqs.csv")

df_mean = df %>% summarize(mean = mean(forward), median = median(forward))
p = ggplot(df) + geom_histogram(aes(x = forward)) + ylab('samples') + xlab('reads') +
  geom_vline(data = df_mean, aes(xintercept = mean), color = "red", linetype = "dashed") +
  theme(axis.text.x = element_text(size = 6)) + scale_x_log10()
ggsave(paste0("/home/mac9jc/cparvum_genomes/plots/read_counts_processed_"
              , format(Sys.time(), "%Y_%m_%d")
              , ".png"), p, width = 8, height = 2, dpi = 600)

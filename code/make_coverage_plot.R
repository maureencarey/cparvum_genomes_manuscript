sink('/home/mac9jc/cparvum_genomes/R_outfile_make_coverage_plot.txt')
files <- list.files(path="/scratch/mac9jc/genome_analysis_project/results",
                    pattern="*coverage.txt", full.names=TRUE, recursive=FALSE)

x = files[1]
coverage=read.table(x, sep="\t", header=F)
colnames(coverage)[colnames(coverage) == 'V1'] = "Chr"
colnames(coverage)[colnames(coverage) == 'V2'] = "locus"
colnames(coverage)[colnames(coverage) == 'V3'] = "depth"
x_short = substring(x, 49,59)
coverage$sample = x_short 
cov_all = coverage

for (x in files[c(-1)]) {
  print(x)
  coverage=read.table(x, sep="\t", header=F)
  colnames(coverage)[colnames(coverage) == 'V1'] = "Chr"
  colnames(coverage)[colnames(coverage) == 'V2'] = "locus"
  colnames(coverage)[colnames(coverage) == 'V3'] = "depth"
  x_short = substring(x, 49,59)
  coverage$sample = x_short 
  cov_all = rbind(cov_all,coverage)
}

cov_all$sample = ifelse(substring(cov_all$sample, 1, 1) == 'S',
                        substring(cov_all$sample, 1, 10),
                        cov_all$sample)


library(ggplot2)
library(dplyr)
cov_all_plot1 = cov_all
cov_all_plot1$sample = as.factor(cov_all_plot1$sample)
# replace sample ids
# levels(cov_all_plot1$sample)[levels(cov_all_plot1$sample)==""] = ""

ggplot(cov_all_plot1) + geom_boxplot(aes(x = sample, y = depth), outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) + ylim(0,800) # + scale_y_sqrt()
ggsave("/home/mac9jc/cparvum_genomes/plots/SF1.pdf", height = 4, width = 6.4)

library(dplyr)
# list of sample ids
# use_samples = c()
figure_df = cov_all[cov_all$sample %in% use_samples,]
figure_df$sample = as.factor(figure_df$sample)
# replace sample ids
# levels(figure_df$sample)[levels(figure_df$sample)==""] = ""
print('mean coverage for Bangladesh samples:')
print(mean(figure_df$depth))
UK = cov_all[!(cov_all$sample %in% use_samples),]
print('mean coverage for UK samples:')
print(mean(UK$depth))
ggplot(figure_df) + geom_boxplot(aes(x = sample, y = depth), outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) + scale_y_sqrt(limits = c(0,450))
ggsave("/home/mac9jc/cparvum_genomes/plots/fig1A.pdf", height = 3, width = 3)

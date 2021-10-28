library(ggplot2)
library(dplyr)
library(stringr)
library(zoo)
print('loaded packages')


#sink('/home/mac9jc/cparvum_genomes/R_outfile_make_LD_plot.txt')
file_list = list.files("/scratch/mac9jc/genome_analysis_project/results/", pattern = "^LD.*\\.ld$")
print(file_list)

setwd("/scratch/mac9jc/genome_analysis_project/results/")
for (file in file_list){
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.table(file, header=TRUE)}
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.table(file, header=TRUE)
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)  }
}

# get rolling averages
df = dataset %>% mutate(distance = BP_B- BP_A) %>%
  arrange(distance)
df_temp = df %>% arrange(distance) 
R2_means = zoo::rollapply(df_temp$R2, 1000,mean,by=200)
dist_means = zoo::rollapply(df_temp$distance, 1000,mean,by=200)
df_use = data.frame(dist_means, R2_means)

#plot
p = ggplot(df_use, aes(x = dist_means, y = R2_means)) + 
  geom_point(alpha = .3) + 
  geom_line() +
  xlab('distance between SNPs (bp)')+
  ylab('mean R squared') 
ggsave("/home/mac9jc/cparvum_genomes/plots/LD_plot.pdf",p, width = 6, height = 4, units = "in", dpi = 600)

print('saved figure 1')
cp = df_use

file_list = list.files("/scratch/mac9jc/genome_analysis_project/chominis/", pattern = "^LD.*\\.ld$")
print(file_list)

i = 0
rm(dataset)
setwd("/scratch/mac9jc/genome_analysis_project/chominis/")
for (file in file_list){
  i = i + 1
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.table(file, header=TRUE)
    if (grepl('_0',file)) {
      dataset = cbind(dataset, subset='full')}
    else {
      dataset = cbind(dataset, subset=i)}}
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.table(file, header=TRUE)
    if (grepl('_0',file)) {
      temp_dataset = cbind(temp_dataset, subset='full')}
    else {
      temp_dataset = cbind(temp_dataset, subset=i)}
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)  }
}

print('subsets')
print(unique(dataset$subset))
print('loaded other data')
print(head(dataset))

# get rolling averages
df = dataset %>% mutate(distance = BP_B- BP_A) %>% arrange(distance)
subset_full = df %>% filter(subset != 'full') %>% select(subset) %>% unique()
df_keep = data.frame(dist_means = numeric(), R2_means = numeric(), subset = character())
print('beginning of rolling avgs')
print(head(df))
for (i in 1:length(unique(df$subset))) {
  subset_temp = unique(df$subset)[i]
  df_temp = df %>% filter(subset == subset_temp) %>% arrange(distance)
  R2_means = zoo::rollapply(df_temp$R2, 1000,mean,by=200)
  dist_means = zoo::rollapply(df_temp$distance, 1000,mean,by=200)
  df_use = data.frame(dist_means, R2_means) %>% mutate(subset = subset_temp)
  df_keep = rbind(df_keep, df_use)
}
print('calc rolling avgs')
print(head(df_keep))
print(unique(df_keep$subset))

# color code full dataset
df_keep = df_keep[-1,] %>% mutate(full_dataset = as.factor(ifelse(subset == "full",1,0)))

p = ggplot(df_keep %>% mutate(subset_use = ifelse(subset == "full", NA, subset)), 
           aes(x = dist_means, y = R2_means, color = subset_use)) +
  geom_point(alpha = .3) +
  geom_line() +
  scale_color_grey(na.value = "blue") +
  xlab('distance between SNPs (bp)')+
  ylab('mean R squared') + guides(color = F)
ggsave("/home/mac9jc/cparvum_genomes/plots/LD_plot_ch_only.pdf",p, width = 6, height = 4, units = "in", dpi = 600)

print('saved figure 2')

ch = df_keep %>% filter(subset == "full") %>% select(-subset,-full_dataset)
print('dimensions')
print(dim(ch))
print(dim(cp))

df = rbind(ch %>% mutate(species = "C. hominis"),
           cp %>% mutate(species = "C. parvum"))
p = ggplot(df, aes(x = dist_means, y = R2_means, color = species)) + 
  geom_point(alpha = .3) + 
  geom_line() +
  xlab('distance between SNPs (bp)')+
  scale_color_manual(values = c('black','blue')) +
  ylab('mean R squared') + xlim(0,5000) +
  geom_vline(xintercept = 300)
ggsave("/home/mac9jc/cparvum_genomes/plots/LD_plot_both_smoother.pdf",p, width = 6, height = 4, units = "in", dpi = 600)

library(tidyverse)
library(vcfR)

files <- list.files(path="~/Documents/Documents\ -\ somH12CXJ8GJV40/work/projects/cparvum_genomes/copy_to_mac/renamed", 
                    full.names=TRUE, recursive=FALSE)
list_o_muts = list()
errors = list()
good= list()
for (f in files) {
  vcf = read.vcfR(f) # load file
  sample = strsplit(basename(f),"\\.")[[1]][2]
  # apply function
  gt = as.data.frame(extract.gt(vcf, element = 'GT'))
  if (("0/0/0/0" %in% gt[,1])) { errors[[sample]] = sample}
  if ((gt[1,] %in% gt[,1])) { good[[sample]] = sample}
  list_o_muts[[sample]] = rownames(gt)
}

# merge all into sparse matrix
prev_df = data.frame(V2=character(),
                     stringsAsFactors=FALSE)
for (x in 1:length(list_o_muts)) {
  nm = names(list_o_muts)[x]
  df = as.data.frame(list_o_muts[nm])
  df[,2] = df[,1]
  prev_df = merge(prev_df, df, by = "V2", all = T)
}

# convert to sparse
rownames(prev_df) = prev_df$V2; prev_df$V2 = NULL
prev_df[!is.na(prev_df)] <- 1
prev_df[is.na(prev_df)] <- 0

# quantify snps and break into groups
all_df = prev_df
colnames(all_df) = gsub(pattern= "_renamed", replacement = "", x = colnames(all_df))
all_df[] <- lapply(all_df, function(x) as.numeric(as.character(x)))

B_df = all_df[,-c(7:18)]; B_df = B_df[apply(B_df, 1, function(x) !all(x==0)),]
UK_df = all_df[,c(7:18)]; UK_df = UK_df[apply(UK_df, 1, function(x) !all(x==0)),]
colSums (all_df, na.rm = FALSE, dims = 1)
dim(B_df)
total_snps_B = dim(B_df)[[1]]
dim(UK_df)
total_snps_UK = dim(UK_df)[[1]]

# visualize
library(UpSetR)
genomes = as.data.frame(readxl::read_xlsx("~/Documents/Documents - somH12CXJ8GJV40/work/projects/cparvum_genomes/genomes.xlsx")) %>%
  mutate(gp60_long = gp60,
         gp60 = ifelse(grepl('mixed infection',gp60_long),'IIc',substr(gp60_long,1,3)),
         type = ifelse(grepl('mixed infection',gp60_long),'mixed infection','monoinfection') ) %>%
  dplyr::rename(site = `isolated in`, sets = `genome name`) %>%
  select(sets, gp60_long,gp60,identifier, type, site) %>%
  mutate(sets = ifelse(substr(sets,1,1) =='U',identifier,sets))

genomes_B = genomes[grepl("Bangladesh",genomes$site),]
genomes_B = genomes_B[,colnames(genomes_B) %in% c("sets", "gp60")]
genomes_B = genomes_B[genomes_B$sets != "icddrb29",]
genomes_UK = genomes[!grepl("Bangladesh",genomes$site),]
genomes_UK = genomes_UK[,colnames(genomes_UK) %in% c("sets", "gp60")]
genomes_B = genomes_B[c("sets","gp60")]
genomes_UK = genomes_UK[c("sets","gp60")]

upset(B_df[,colnames(B_df)!= "icddrb29"], main.bar.color = "black")
pdf(file="~/Documents/Documents - somH12CXJ8GJV40/work/projects/cparvum_genomes/plots/UpsetR_B.pdf", onefile=FALSE, width=5, height=5)
upset(B_df[,colnames(B_df)!= "icddrb29"], 
      sets = c('icddrb63', 'icddrb47', 'icddrb111','icddrb90','icddrb93'),
      intersections = list(
        list("icddrb63"),
        list("icddrb47"),
        list("icddrb111"),
        list("icddrb90"),
        list("icddrb93"),
        list('icddrb63', 'icddrb47', 'icddrb111','icddrb90','icddrb93'), # all plus largest 3, gp60
        list('icddrb47', 'icddrb93'),
        list('icddrb63', 'icddrb47'),
        list('icddrb111','icddrb90','icddrb93'),
        list('icddrb47', 'icddrb111','icddrb90','icddrb93')),
      mainbar.y.label = "number of shared SNP loci",
      sets.x.label = "SNPs",
      set.metadata = list(data = genomes_B, 
                          plots = list(list(type = "matrix_rows", column = "gp60", colors = c(IIc = "#E69F00", IId = "#0072B2"), alpha = 0.5))),
      keep.order = TRUE)
dev.off()

upset(UK_df, main.bar.color = "black", nsets = 12)
pdf(file="~/Documents/Documents - somH12CXJ8GJV40/work/projects/cparvum_genomes/plots/UpsetR_UK.pdf", onefile=FALSE, width=5, height=5)
upset(UK_df, main.bar.color = "black", 
      intersections = list(list("SRR6117460", "SRR6147472", "SRR6147581", "SRR6147587", "SRR6147945",
                                "SRR6147964", "SRR6148259", "SRR6813716", "SRR6813717", "SRR6813718",
                                "SRR6813719", "SRR7898459"),
                           list("SRR6117460"),
                           list("SRR6147472"),
                           list("SRR6147581"),
                           list("SRR6147587"),
                           list("SRR6147945"),
                           list("SRR6147964"),
                           list("SRR6148259"),
                           list("SRR6813716"),
                           list("SRR6813717"),
                           list("SRR6813718"),
                           list("SRR6813719"),
                           list("SRR7898459"),
                           list("SRR6813716",  "SRR6813718","SRR6813719"),
                           list("SRR6813718","SRR6813719"),
                           list("SRR6813716",  "SRR6813718"),
                           list("SRR6148259", "SRR7898459"),
                           list("SRR6117460", "SRR6147472", "SRR6147581", "SRR6147587", "SRR6147945",
                                "SRR6147964"),
                           list("SRR6813716", "SRR6813717", "SRR6813718",
                                "SRR6813719", "SRR7898459")),
      mainbar.y.label = "number of shared SNP loci",
      sets.x.label = "SNPs",
      set.metadata = list(data = genomes_UK, 
                          plots = list(list(type = "matrix_rows", column = "gp60", colors = c(IIc = "#E69F00", IId = "#0072B2", IIa = "#000000"), alpha = 0.5))))
dev.off()

# get snps in all of the B samples/ UK samples
all_B = B_df[,colnames(B_df)!= "icddrb29"]
all_B = all_B[rowSums(all_B)==ncol(all_B),]
all_UK = UK_df[rowSums(UK_df)==ncol(UK_df),]
df = merge(all_UK %>% select("SRR6117460"), all_B %>% select("icddrb47"), by = "row.names", all = T)
df[is.na(df)] <- 0
colnames(df) = c("loci","UK","Bangladesh")
pdf(file="~/Documents/Documents - somH12CXJ8GJV40/work/projects/cparvum_genomes/plots/UpsetR_B_UK.pdf", onefile=FALSE, width=3, height=5)
upset(df[,-1], main.bar.color = "black",
      mainbar.y.label = "number of shared SNP loci",
      sets.x.label = "SNPs")
dev.off()


###### SNP LOCATION ANALYSIS

files <- list.files(path="~/Documents/Documents\ -\ somH12CXJ8GJV40/work/projects/cparvum_genomes/copy_to_mac/all", 
                    pattern="*_2.vcf", full.names=TRUE, recursive=FALSE)
list_o_muts = list()
errors = list()
good= list()
for (f in files) {
  vcf = read.vcfR(f) # load file
  sample = strsplit(basename(f),"\\.")[[1]][2]
  # apply function
  gt = as.data.frame(extract.gt(vcf, element = 'GT'))
  if (("0/0/0/0" %in% gt[,1])) { errors[[sample]] = sample}
  if ((gt[1,] %in% gt[,1])) { good[[sample]] = sample}
  list_o_muts[[sample]] = rownames(gt)
}

# merge all into sparse matrix
prev_df = data.frame(V2=character(),
                     stringsAsFactors=FALSE)
for (x in 1:length(list_o_muts)) {
  nm = names(list_o_muts)[x]
  df = as.data.frame(list_o_muts[nm])
  df[,2] = df[,1]
  prev_df = merge(prev_df, df, by = "V2", all = T)
}

# convert to sparse
rownames(prev_df) = prev_df$V2; prev_df$V2 = NULL
prev_df[!is.na(prev_df)] <- 1
prev_df[is.na(prev_df)] <- 0

# quantify snps and break into groups
all_df = prev_df
all_df[] <- lapply(all_df, function(x) as.numeric(as.character(x)))
B_df = all_df[,-c(7:18)]; B_df = B_df[apply(B_df, 1, function(x) !all(x==0)),]
UK_df = all_df[,c(7:18)]; UK_df = UK_df[apply(UK_df, 1, function(x) !all(x==0)),]
colSums (all_df, na.rm = FALSE, dims = 1)
dim(B_df)
# percent missense
dim(B_df)/total_snps_B
dim(UK_df)
# percent missense
dim(UK_df)/total_snps_UK

# get snps in all of the B samples/ UK samples MISSENSE ONLY
all_B = B_df[,colnames(B_df)!= "icddrb29_renamed"]
all_B = all_B[rowSums(all_B)==ncol(all_B),]
all_UK = UK_df[rowSums(UK_df)==ncol(UK_df),]
df = merge(all_UK %>% select("SRR6117460_renamed"), all_B %>% select("icddrb47_renamed"), by = "row.names", all = T)
df[is.na(df)] <- 0
colnames(df) = c("loci","UK","Bangladesh")
pdf(file="~/Documents/Documents - somH12CXJ8GJV40/work/projects/cparvum_genomes/plots/UpsetR_B_missense.pdf", onefile=FALSE, width=5, height=5)
upset(df[,-1], main.bar.color = "black",
      mainbar.y.label = "number of shared SNP loci",
      sets.x.label = "SNPs")
dev.off()

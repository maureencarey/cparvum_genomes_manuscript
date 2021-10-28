
# make sure VCF and metadata files are here
setwd("~/Documents/Documents - somH12CXJ8GJV40/work/projects/cparvum_genomes/")

# load packages that contain the functions used below
library(tidyverse)
library(readxl)
SEED_INT = 2019
set.seed(SEED_INT)


##### ALL SNPS #####

pca <- read_table("~/Documents/Documents - somH12CXJ8GJV40/work/projects/cparvum_genomes/pca_results/pca.eigenvec", col_names = FALSE)
eigenval <- scan("~/Documents/Documents - somH12CXJ8GJV40/work/projects/cparvum_genomes/pca_results/pca.eigenval")

# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
pca_all = pca

# get relevant groupings
location <- rep(NA, length(pca$ind))
location[grep("SRR", pca$ind)] <- "UK"
location[grep("__", pca$ind)] <- "Bangladesh"
pca <- as_tibble(data.frame(pca, location)) 
genomes = readxl::read_xlsx("genomes.xlsx") %>%
  mutate(gp60_long = gp60,
         gp60 = ifelse(grepl('mixed infection',gp60_long),'IIc',substr(gp60_long,1,3)),
         type = ifelse(grepl('mixed infection',gp60_long),'mixed infection','monoinfection') ) %>%
  select(gp60_long,gp60,identifier, type)
pca = merge(pca, genomes, by.x = "ind", by.y = "identifier", all = T)
rm(location)

# first convert to percentage variance explained
pve = data.frame(PC = 1:18, pve = eigenval/sum(eigenval)*100)
a = ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# plot pca
ggplot(pca, aes(PC1, PC2, col = gp60, shape = location, size = type)) + 
  geom_point(position = position_jitter(width = 0.02,height = 0.02), alpha = 0.6) +
  coord_equal() + theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + 
  scale_color_manual(values = c("darkgrey","#E69F00","#0072B2")) #+
  #ggrepel::geom_label_repel(aes(label = ind))
ggsave("~/Documents/Documents - somH12CXJ8GJV40/work/projects/cparvum_genomes/plots/pca_all.png",
       width = 9,height = 9,units = "in",dpi = 300)

vegan::adonis(pca[,2:19] ~ pca$location + pca$gp60, method="euclidean",perm=999)

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# pca$location  1         1       1       1 0.05882  0.550  
# pca$gp60      2         2       1       1 0.11765  0.034 *
#   Residuals    14        14       1         0.82353         
# Total        17        17                 1.00000   

##### MISSENSE ONLY #####
dev.off()

pca <- read_table("~/Documents/Documents - somH12CXJ8GJV40/work/projects/cparvum_genomes/pca_results/pca_missense.eigenvec", col_names = FALSE)
eigenval <- scan("~/Documents/Documents - somH12CXJ8GJV40/work/projects/cparvum_genomes/pca_results/pca_missense.eigenval")

# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# get relevant groupings
location <- rep(NA, length(pca$ind))
location[grep("SRR", pca$ind)] <- "UK"
location[grep("__", pca$ind)] <- "Bangladesh"
pca <- as_tibble(data.frame(pca, location)) 
pca = merge(pca, genomes, by.x = "ind", by.y = "identifier", all = T)
rm(location)

# plot pca
fig_A_leg = ggplot(pca, aes(PC1, PC2, col = gp60, shape = location, size = type)) + 
  geom_point(position = position_jitter(width = 0.02,height = 0.02), alpha = 0.5) +
  coord_equal() + theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + 
  scale_color_manual(values = c("darkgrey","#E69F00","#0072B2")) +
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         shape=guide_legend(nrow=1, byrow=TRUE),
         size=guide_legend(nrow=2, byrow=TRUE)) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))#+
#ggrepel::geom_label_repel(aes(label = ind))
fig_A = ggplot(pca, aes(PC1, PC2, col = gp60, shape = location, size = type)) + 
  geom_point(position = position_jitter(width = 0.02,height = 0.02), alpha = 0.5) +
  coord_equal() + theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + 
  scale_color_manual(values = c("darkgrey","#E69F00","#0072B2")) + 
  theme(legend.position="none", plot.margin=unit(c(1,1,1,1),"cm"))
fig_A


vegan::adonis(pca[,2:19] ~ pca$location + pca$gp60, method="euclidean",perm=999)
# NEITHER ARE SIGNIFICANT


##### MISSENSE SECRETED ONLY #####
dev.off()

pca <- read_table("~/Documents/Documents - somH12CXJ8GJV40/work/projects/cparvum_genomes/pca_results/pca_missense_secreted.eigenvec", col_names = FALSE)
eigenval <- scan("~/Documents/Documents - somH12CXJ8GJV40/work/projects/cparvum_genomes/pca_results/pca_missense_secreted.eigenval")


##### MISSENSE SECRETED ONLY WITHOUT GP60 #####

dev.off()

pca <- read_table("~/Documents/Documents - somH12CXJ8GJV40/work/projects/cparvum_genomes/pca_results/pca_missense_secreted_no_gp60.eigenvec", col_names = FALSE)
eigenval <- scan("~/Documents/Documents - somH12CXJ8GJV40/work/projects/cparvum_genomes/pca_results/pca_missense_secreted_no_gp60.eigenval")

# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# get relevant groupings
location <- rep(NA, length(pca$ind))
location[grep("SRR", pca$ind)] <- "UK"
location[grep("__", pca$ind)] <- "Bangladesh"
pca <- as_tibble(data.frame(pca, location)) 
pca = merge(pca, genomes, by.x = "ind", by.y = "identifier", all = T)
rm(location)

# plot pca
fig_B = ggplot(pca, aes(PC1, PC2, col = gp60, shape = location, size = type)) + 
  geom_point(position = position_jitter(width = 0.02,height = 0.02), alpha = 0.5) +
  coord_equal() + theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + 
  scale_color_manual(values = c("darkgrey","#E69F00","#0072B2"))+ 
  theme(legend.position="none", plot.margin=unit(c(1,1,1,1),"cm"))#+
#ggrepel::geom_label_repel(aes(label = ind)) # OUTLIER IS 47
fig_B

vegan::adonis(pca[,2:19] ~ pca$location + pca$gp60, method="euclidean",perm=999)
# neither are significant


##### NO OUTLIER #####SRR6813716 didn't pass mind

dev.off()

pca <- read_table("~/Documents/Documents - somH12CXJ8GJV40/work/projects/cparvum_genomes/pca_results/no_47_pca.eigenvec", col_names = FALSE)
eigenval <- scan("~/Documents/Documents - somH12CXJ8GJV40/work/projects/cparvum_genomes/pca_results/no_47_pca.eigenval")

# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# get relevant groupings
location <- rep(NA, length(pca$ind))
location[grep("SRR", pca$ind)] <- "UK"
location[grep("__", pca$ind)] <- "Bangladesh"
pca <- as_tibble(data.frame(pca, location)) 
# remove_sample = # set this to sample id
pca = merge(pca, genomes %>% filter(identifier != remove_sample, identifier != "SRR6813716"), by.x = "ind", by.y = "identifier", all = T)
rm(location)

# plot pca
fig_C = ggplot(pca, aes(PC1, PC2, col = gp60, shape = location, size = type)) + 
  geom_point(position = position_jitter(width = 0.02,height = 0.02), alpha = 0.5) +
  coord_equal() + theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + 
  scale_color_manual(values = c("darkgrey","#E69F00","#0072B2")) + 
  theme(legend.position="none", plot.margin=unit(c(1,1,1,1),"cm"))#+
#ggrepel::geom_label_repel(aes(label = ind))
fig_C

vegan::adonis(pca[,2:17] ~ pca$location + pca$gp60, method="euclidean",perm=999)
# neither are significant

leg = ggpubr::as_ggplot(ggpubr::get_legend(fig_A_leg))
gridExtra::grid.arrange(fig_A,fig_B,fig_C,leg, nrow = 2)
#ggpubr::ggarrange(fig_A,fig_B,fig_C,ncol=3, nrow=1, common.legend = TRUE, legend="bottom")

ggsave("~/Documents/Documents - somH12CXJ8GJV40/work/projects/cparvum_genomes/plots/extra_PCAs.png",
       width = 9,height = 9,units = "in",dpi = 300)

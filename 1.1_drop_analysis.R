library(dartRverse)
library(dartR.popgen)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
source("scripts/theme_emily.R")

#~~ Script to run downsampled data

# Write metadata file for subsampled individuals

fread("data/out/plink/gl_sec_het_rhina_drop.fam") %>%
  select(V2, V1) %>%
  rename(id = V2, pop = V1) %>%
  write.csv("data/meta/rhina_sub_metadata.csv", quote = F, 
            row.names = F)

metadata <- "data/meta/rhina_sub_metadata.csv"

gl_drop <- gl.read.PLINK(
  "data/out/plink/gl_sec_het_rhina_drop",
  ind.metafile = metadata,
  loc.metafile = NULL,
  plink.cmd = "plink")

gl_drop <- gl.reassign.pop(gl_drop, as.pop="Family", verbose = NULL)

nLoc(gl_drop) # 3463
nInd(gl_drop) # 24
gl_drop@pop

#~~~~~~~~~~~~~~~~~~~~~~#
#         PCA          #
#~~~~~~~~~~~~~~~~~~~~~~#

gl_drop_maf <- gl.filter.maf(gl_drop, threshold=0.03, verbose=3) # equivalent to MAC 3

nLoc(gl_drop_maf) # 3193
nInd(gl_drop_maf) # 24

pca1 <- glPca(gl_drop_maf, nf = nLoc(gl_drop_maf))

pc1 <- pca1$scores[,1]
pc2 <- pca1$scores[,2]
pc3 <- pca1$scores[,3]
pc4 <- pca1$scores[,4]

ind_names <- gl_drop_maf@ind.names
pop_names <- gl_drop_maf@pop

ggplot_pca <- as.data.frame(cbind(pc1, pc2, pc3, pc4)) %>%
  mutate(ind_names = ind_names) %>%
  mutate(pop = pop_names) %>%
  mutate(pop = case_when(grepl("Sri_Lanka", pop) ~ "Sri Lanka",
                         grepl("Saudi_Arabia", pop) ~ "Saudi Arabia",
                         TRUE ~ pop))
  

eig <- data.frame(pca1$eig)
eig$percentage = (eig[, 1]/sum(eig$pca1.eig))*100
sum(eig$percentage)
sum(eig$percentage[1:2])

eig$percentage <- round(eig$percentage, digits = 1)
eig$percentage[1]
eig$percentage[2]
eig$percentage[3]

col_palette <- c("#F2B705",
                 "#046C9A",
                 "#FD6467",
                 "#9986A5",
                 "#0B775E",
                 "#D67236",
                 "#4F1C2E")


pc1_pc2 <- ggplot(ggplot_pca, aes(pc1, pc2)) + 
  geom_point(aes(col = factor(pop)), size = 4, alpha = 0.7) +
  theme_emily() +
  scale_color_manual(values = col_palette, name = "Location") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain"),
        legend.position = "none") +
  ylab(paste0("PC2 (",eig$percentage[2],"%)")) +
  xlab(paste0("PC1 (",eig$percentage[1],"%)")) + 
  labs(title = "A")

pc1_pc3 <- ggplot(ggplot_pca, aes(pc1, pc3)) + 
  geom_point(aes(col = factor(pop)), size = 4, alpha = 0.7) +
  theme_emily() +
  scale_color_manual(values = col_palette, name = "Location") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain"),
        plot.title = element_text(face = "plain")) +
  ylab(paste0("PC3 (",eig$percentage[3],"%)")) +
  xlab(paste0("PC1 (",eig$percentage[1],"%)")) +
  labs(title = "B")
  

pc1_pc2 + pc1_pc3


ggsave("figs/PCA_rhina_sub.png", pc1_pc2 / pc1_pc3, height = 7, width = 6)
ggsave("figs/PCA_rhina_sub_12.png", pc1_pc2, height = 3, width = 5)
ggsave("figs/PCA_rhina_sub_13.png", pc1_pc3, height = 4, width = 6)


saveRDS(pc1_pc2, "figs/PCA_rhina_sub_12.RDS")


k2_plot <- readRDS("figs/structure_sub_K2_plot.RDS")
k3_plot <- readRDS("figs/structure_sub_K3_plot.RDS")
k4_plot <- readRDS("figs/structure_sub_K4_plot.RDS")
k5_plot <- readRDS("figs/structure_sub_K5_plot.RDS")
k6_plot <- readRDS("figs/structure_sub_K6_plot.RDS")
k7_plot <- readRDS("figs/structure_sub_K7_plot.RDS")
k8_plot <- readRDS("figs/structure_sub_K8_plot.RDS")

deltak <- readRDS("figs/deltak_sub.RDS")
loglike <- readRDS("figs/loglike_sub.RDS")


fig_s7 <- ((deltak / loglike + plot_layout(guides = 'auto')) | 
             k2_plot / k3_plot / k4_plot / k5_plot / k6_plot) + 
  plot_layout(guides = 'collect', widths = c(1.5, 1))

fig_s7

ggsave("figs/fig_s7.png", fig_s7, height = 7, width = 7)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#            Fst               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

gl_drop_maf@pop <- as.factor(ggplot_pca$pop)

fst <- gl.fst.pop(gl_drop_maf, nboots=1000, percent=95, nclusters=1)
pw_fst_boot <- fst$Bootstraps[,c(1,2,1003:1006)]

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}


fst$Fsts[is.na(fst$Fsts)] <- t(fst$Fsts)[is.na(fst$Fsts)]

fst$Pvalues

# melted_cor <- fst$Fsts
melted_cor <- fst$Fsts[c(1,2,3),c(1,2,3)] # Remove saudi and bangladesh due to small sample sizes
melted_cor <- fst$Fsts[c(4,2,1,3,5),c(4,2,1,3,5)] # Reorder geographically

upper_tri <- get_upper_tri(melted_cor)
melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)

fst_plot <- ggplot(melted_cormat, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  #scale_fill_gradient2(name="Fst") +
  scale_fill_gradient(low = "#F8FCFC",
                      high = "#35274A",
                      name = expression(italic("F")[ST])) + #limit = c(-0.0025,0.160)
  theme_emily() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  coord_fixed()

fst_plot


ggsave("figs/FST_rhina_sub_SNPs.png", fst_plot, height = 4, width = 5)


# Plotting CIs

fst_cis <- pw_fst_boot %>%
  unite(Pairwise, c("Population1", "Population2")) %>%
  mutate(Pairwise = fct_reorder(Pairwise, desc(Fst))) %>%
  mutate(sig = as.factor(case_when(`Lower bound CI limit` < 0 & `Upper bound CI limit` > 0 ~ 1,
                                   TRUE ~ 0))) %>%
  filter(!grepl("AP", Pairwise)) %>%
  mutate(col =  "black") %>%
  mutate(model = paste0("<span style=\"color: ", col, "\">", Pairwise, "</span>")) %>%
  ggplot(aes(Fst, fct_reorder(model, Fst), col = sig, fill = sig)) +
  geom_pointrange(aes(xmax = `Upper bound CI limit`, xmin = `Lower bound CI limit`), 
                  size = 0.3, 
                  position=position_dodge(0.2)) +
  geom_vline(xintercept = 0) +
  scale_color_manual(values = c("black", "grey70"), name = "Significant") +
  theme_emily() +
  theme(legend.position = "none",
        plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12),
        axis.text.y = ggtext::element_markdown()) +
  ggtitle("B") + xlab(expression(italic("F")[ST]))

fst_cis

ggsave("figs/FST_panel_rhina_drop.png", fst_plot + fst_cis, height = 6, width = 10)


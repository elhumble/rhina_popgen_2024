library(dartR)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
source("scripts/theme_emily.R")

dartfile <- "data/raw/OrderAppendix_2_DRhin22-7571/Report_DRhin22-7571_SNP_1.csv"
metadata <- "data/meta/rhina_metadata.csv"

gl <- gl.read.dart(dartfile, ind.metafile = metadata, probar=TRUE)

# drop 4 additional species

gl <- gl.drop.ind(
  gl,
  ind.list = c("ID25", 
               "ID67", 
               "ID68", 
               "ID69"))


# metrics

hist(gl@other$loc.metrics$CallRate)
nLoc(gl) # 19436
nInd(gl) # 65
gl.smearplot(gl)

#~~ filter on reproducibility 

gl.report.reproducibility(gl) 

gl <- gl.filter.reproducibility(gl,
                                threshold = 0.99, verbose = 3)

# 3060 loci discarded, 16376 loci remaining

gl <- gl.recalc.metrics(gl)

#~~ select one SNP per sequence tag

gl_sec <- gl.filter.secondaries(gl, verbose = 3, method = "best")

# 892 SNPs removed, 15484 remaining 

#~~ depth

gl.report.rdepth(gl_sec)

# lower 5 and 50

gl_sec <- gl.filter.rdepth(gl_sec, lower = 5,
                           upper = 50)

nLoc(gl_sec) # 13062 remaining

#~~ call rate

# SNP call rate
gl.report.callrate(gl_sec)

# ind call rate
gl.report.callrate(gl_sec, method="ind")

# filter snp call rate

gl_sec <- gl.filter.callrate(gl_sec, method="loc", 
                             threshold = 0.80,
                             verbose = 3)

# 12576 / 13062 loci retained

gl.report.callrate(gl_sec, method="ind")

# filter ind call rate

gl_sec <- gl.filter.callrate(gl_sec, method="ind", 
                             threshold = 0.96,
                             verbose = 3)

nLoc(gl_sec)
nInd(gl_sec)

# 65 / 69 individuals retained
# ID67 / ID68 / ID69 removed (all bangladesh) / ID35 oman

gl.smearplot(gl_sec)

# Individual heterozygosity

ind_het <- gl.report.heterozygosity(
  gl_sec,
  method = "ind",
  verbose = 3
)

# Low vol samples
# ID08
# ID20
# ID28
# ID32
# ID29

# High vol samples
# ID48

het_out <- arrange(filter(ind_het, Ho > 0.09), Ho)


# filter for major outlier heterozygosity and high and low vol samples

gl_sec_het <- gl.drop.ind(gl_sec,
                          ind.list=c("ID28",
                                     "ID20",
                                     "ID22",
                                     "ID24",
                                     "ID29",
                                     "ID58",
                                     "ID19",
                                     "ID62",
                                     "ID54",
                                     "ID08",
                                     "ID32",
                                     "ID48"
                          ))

nLoc(gl_sec_het) # 12576
nInd(gl_sec_het) # 52

gl_sec_het <- gl.filter.monomorphs(gl_sec_het)
gl_sec_het <- gl.recalc.metrics(gl_sec_het) 

nLoc(gl_sec_het) # 6391
nInd(gl_sec_het) # 52

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Write PLINK files for relatedness      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Write plink files for relatedness

chr <- tibble(locnames = gl_sec_het@loc.names) %>%
  separate(locnames, c("chr", "pos", "geno"), sep = "-") %>%
  mutate(chr = as.factor(chr),
         pos = as.numeric(pos))

gl_plink <- gl_sec_het

gl_plink@chromosome <- chr$chr

dartR::gl2plink(gl_plink,
                plink_path = "/Users/emilyhumble/software",
                bed_file = T,
                outfile = "relatedness",
                outpath = "data/out/plink/",
                pos_cM = "0",
                ID_dad = "0",
                ID_mom = "0",
                sex_code = "unknown",
                phen_value = "0",
                verbose = 3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#         Relatedness         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

infile <- "data/out/plink/relatedness"
system(paste0("/Users/emilyhumble/software/plink --bfile ", infile, " --genome --maf 0.1 --geno 0.1 --hwe 0.001 --recode vcf-iid --allow-extra-chr --debug --out data/out/plink/relatedness"))

# Run NGSrelate
system(paste0("/Users/emilyhumble/software/ngsRelate/ngsRelate -h data/out/plink/relatedness.vcf -T GT -O data/out/plink/relatedness.res -c 1"))

# Read in output

gen <- fread("data/out/plink/relatedness.genome", header = T)
ngsrel <- fread("data/out/plink/relatedness.res")

gen$R1 <- ngsrel$R1 
gen$R0 <- ngsrel$R0
gen$KING <- ngsrel$KING

gen <- gen %>%
  mutate(kinship = (Z1 / 4) + (Z2 / 2)) %>%
  mutate(criteria = case_when(kinship >= 1/2^(5/2) & kinship < 1/2^(3/2) & Z0 < 0.1 ~ "Parent-offspring",
                              kinship >= 1/2^(5/2) & kinship < 1/2^(3/2) & Z0 > 0.1 & Z0 < 0.365 ~ "Full-sibling",
                              kinship >= 1/2^(7/2) & kinship < 1/2^(5/2) & Z0 > 0.365 & Z0 < 1-(1/(2^(3/2))) ~ "Second-degree",
                              kinship >= 1/2^(9/2) & kinship < 1/2^(7/2) & Z0 > 1-(1/2^(3/2)) & Z0 < 1 -(1/2^(5/2)) ~ "Third-degree",
                              kinship < 1/2^(9/2) & Z0 > 1-(1/2^(5/2)) ~ "Unrelated",
                              TRUE ~ "Unknown"))

summary(gen$PI_HAT)
hist(gen$PI_HAT)

relate_fig <- ggplot(gen, aes(R1, KING)) +
  geom_point(size = 2, alpha = 0.3, shape = 21,
             aes(fill = ifelse(KING > 1/2^(5/2), "red", "black"),
                 col = ifelse(KING > 1/2^(5/2), "red", "black"))) +
  scale_colour_manual(values=c("black", "red")) + 
  scale_fill_manual(values=c("black", "red")) + 
  theme_emily() +
  scale_x_log10() +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(alpha=1),
                               keyheight = 0.1,
                               default.unit = "inch")) +
  geom_hline(yintercept = 1/2^(5/2), linetype = "dashed", alpha = 0.3) #+
  #geom_vline(xintercept = 0.9, linetype = "dashed", alpha = 0.3)

relate_fig

# filter one individual from each of the 3 related pairs

ggsave("figs/relatedness_rhina.png", relate_fig, height = 4, width = 5)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#         Further filtering         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

gl.report.callrate(gl_sec_het, method="ind")


# drop related individuals

gl_sec_het_rel <- gl.drop.ind(gl_sec_het,
                              ind.list=c("ID49",
                                         "ID21",
                                         "ID31"
                              ))

nInd(gl_sec_het_rel)
nLoc(gl_sec_het_rel)

# filter monomorphic loci

gl_sec_het_rel <- gl.recalc.metrics(gl_sec_het_rel) 
gl_sec_het_rel <- gl.filter.monomorphs(gl_sec_het_rel)

nInd(gl_sec_het_rel) # 49
nLoc(gl_sec_het_rel) # 6284

# Check het

ind_het <- gl.report.heterozygosity(
  gl_sec_het_rel,
  method = "ind",
  verbose = 3
)


nLoc(gl_sec_het_rel) # 6284
nInd(gl_sec_het_rel) # 49

#~~~~~~~~~~~~~~~~~~~~~~#
#         PCA          #
#~~~~~~~~~~~~~~~~~~~~~~#

gl_sec_het_rel_maf <- gl.filter.maf(gl_sec_het_rel, threshold=0.03, verbose=3) # equivalent to MAC 3

nLoc(gl_sec_het_rel_maf) # 3533
nInd(gl_sec_het_rel_maf) # 49

pca1 <- glPca(gl_sec_het_rel_maf, nf = nLoc(gl_sec_het_rel_maf))

pc1 <- pca1$scores[,1]
pc2 <- pca1$scores[,2]
pc3 <- pca1$scores[,3]
pc4 <- pca1$scores[,4]

ind_names <- gl_sec_het_rel_maf@ind.names
pop_names <- gl_sec_het_rel_maf@pop

ggplot_pca <- as.data.frame(cbind(pc1, pc2, pc3, pc4)) %>%
  mutate(ind_names = ind_names) %>%
  mutate(pop = pop_names) %>%
  mutate(location = case_when(grepl("oman", pop) ~ "oman",
                              grepl("saudi", pop) ~ "saudi_arabia",
                              grepl("sharjah", pop) ~ "uae",
                              grepl("abudhabi", pop) ~ "uae",
                              grepl("rak", pop) ~ "uae",
                              grepl("dubai", pop) ~ "uae",
                              grepl("srilanka", pop) ~ "srilanka",
                              grepl("bangladesh", pop) ~ "bangladesh"))


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


ggplot_pca <- ggplot_pca %>%
  mutate(location = case_when(location == "uae" ~ "UAE",
                              location == "oman" ~ "Oman",
                              location == "srilanka" ~ "Sri Lanka",
                              location == "saudi_arabia" ~ "Saudi Arabia",
                              location == "bangladesh" ~ "Bangladesh"))



pc1_pc2 <- ggplot(ggplot_pca, aes(pc1, pc2)) + 
  geom_point(aes(col = factor(location)), size = 4, alpha = 0.7) +
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
  geom_point(aes(col = factor(location)), size = 4, alpha = 0.7) +
  theme_emily() +
  scale_color_manual(values = col_palette, name = "Location") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain"),
        plot.title = element_text(face = "italic")) +
  ylab(paste0("PC3 (",eig$percentage[3],"%)")) +
  xlab(paste0("PC1 (",eig$percentage[1],"%)"))

pc1_pc2 + pc1_pc3

ggplot(ggplot_pca, aes(pc2, pc3)) + 
  geom_point(aes(col = factor(location)), size = 4, alpha = 0.7) +
  theme_emily() +
  scale_color_manual(values = col_palette, name = "Location") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain"),
        plot.title = element_text(face = "italic")) +
  ylab(paste0("PC3 (",eig$percentage[3],"%)")) +
  xlab(paste0("PC2 (",eig$percentage[1],"%)"))


ggsave("figs/PCA_rhina.png", pc1_pc2 + pc1_pc3, height = 4, width = 11)
ggsave("figs/PCA_rhina_12.png", pc1_pc2, height = 3, width = 5)
ggsave("figs/PCA_rhina_13.png", pc1_pc3, height = 4, width = 6)


saveRDS(pc1_pc2, "figs/PCA_rhina_12.RDS")

#~~~~~~~~~~~~~~~~~~~~~~~#
#         DAPC          #
#~~~~~~~~~~~~~~~~~~~~~~~#

# Run a DAPC using site IDs as priors

gl_sec_het_rel_maf@pop <- as.factor(ggplot_pca$location)

dapc2 <- dapc(gl_sec_het_rel_maf,
             gl_sec_het_rel_maf$pop,
             n.pca = 4,
             n.da = 3)


percent <- dapc2$eig/sum(dapc2$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,60),
        names.arg = round(percent, 1))

# Create a data.frame containing individual coordinates
ind_coords <- as.data.frame(dapc2$ind.coord)

# Rename columns of dataframe
colnames(ind_coords) <- c("Axis1", "Axis2", "Axis3")

# Add individual names
ind_coords$Ind <- indNames(gl_sec_het_rel_maf)

# Add locations
ind_coords$Site <- gl_sec_het_rel_maf$pop

# Calculate centroid (average) position for each population
centroid <- aggregate(cbind(Axis1, Axis2, Axis3) ~ 
                        Site, data = ind_coords, FUN = mean)

centroid$Site <- factor(centroid$Site)

# Add centroid coordinates to ind_coords dataframe
ind_coords <- left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))

# Custom x and y labels
xlab <- paste("PC1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab <- paste("PC2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Plotting
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2)) +
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), 
               show.legend = FALSE) +
  geom_point(aes(fill = Site, col = Site), size = 3, alpha = 0.7) +
  #geom_label(data = centroid, aes(label = Site, fill = Site, 
  #                                alpha = 0.6), size = 3, show.legend = FALSE) +
  scale_fill_manual(values = col_palette, name = "Location") +
  scale_colour_manual(values = col_palette, name = "Location") +
  labs(x = xlab, y = ylab) +
  #ggtitle("Mobula alfredi") +
  # stat_ellipse(aes(Axis1, Axis2, colour = Site),
  #              type = "norm",
  #              lty = 2, size = 1) +
  theme_emily()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#            Fst               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

gl_sec_het_rel_maf@pop <- as.factor(ggplot_pca$location)

fst <- dartR::gl.fst.pop(gl_sec_het_rel_maf, nboots=1000, percent=95, nclusters=1)
pw_fst_boot_alf <- fst$Bootstraps[,c(1,2,1003:1006)]

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
melted_cormat <- melt(upper_tri, na.rm = TRUE)

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


ggsave("figs/FST_rhina_SNPs.png", fst_plot, height = 4, width = 5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#           Admixture           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# edit gl object prior to converting to plink

chr <- tibble(locnames = gl_sec_het_rel_maf@loc.names) %>%
  separate(locnames, c("chr", "pos", "geno"), sep = "-") %>%
  mutate(chr = as.factor(chr),
         pos = as.numeric(pos))

gl_plink <- gl_sec_het_rel_maf

gl_plink@chromosome <- chr$chr

dartR::gl2plink(gl_plink,
                plink_path = "/Users/emilyhumble/software",
                bed_file = T,
                outfile = "gl_sec_het_rhina",
                outpath = "data/out/plink",
                pos_cM = "0",
                ID_dad = "0",
                ID_mom = "0",
                sex_code = "unknown",
                phen_value = "0",
                verbose = 3)

#~~ Run admixture

system("/Users/emilyhumble/software/admixture_macosx-1.3.0/admixture --cv -B1000 /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/plink/gl_sec_het_rhina.bed 1 -j4 > /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/admixture/gl_sec_het_rhina_admixture_1.out")
system("/Users/emilyhumble/software/admixture_macosx-1.3.0/admixture --cv -B1000 /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/plink/gl_sec_het_rhina.bed 2 -j4 > /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/admixture/gl_sec_het_rhina_admixture_2.out")
system("/Users/emilyhumble/software/admixture_macosx-1.3.0/admixture --cv -B1000 /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/plink/gl_sec_het_rhina.bed 3 -j4 > /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/admixture/gl_sec_het_rhina_admixture_3.out")
system("/Users/emilyhumble/software/admixture_macosx-1.3.0/admixture --cv -B1000 /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/plink/gl_sec_het_rhina.bed 4 -j4 > /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/admixture/gl_sec_het_rhina_admixture_4.out")
system("/Users/emilyhumble/software/admixture_macosx-1.3.0/admixture --cv -B1000 /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/plink/gl_sec_het_rhina.bed 5 -j4 > /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/admixture/gl_sec_het_rhina_admixture_5.out")
system("/Users/emilyhumble/software/admixture_macosx-1.3.0/admixture --cv -B1000 /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/plink/gl_sec_het_rhina.bed 6 -j4 > /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/admixture/gl_sec_het_rhina_admixture_6.out")
system("/Users/emilyhumble/software/admixture_macosx-1.3.0/admixture --cv -B1000 /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/plink/gl_sec_het_rhina.bed 7 -j4 > /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/admixture/gl_sec_het_rhina_admixture_7.out")
system("/Users/emilyhumble/software/admixture_macosx-1.3.0/admixture --cv -B1000 /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/plink/gl_sec_het_rhina.bed 8 -j4 > /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/admixture/gl_sec_het_rhina_admixture_8.out")


# output files always saved in working directory, move to /admixture

system("mv gl_sec* data/out/admixture/")
system("grep 'CV' data/out/admixture/gl*out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' > data/out/admixture/cv_error")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#     Downsample individuals      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

set.seed(123)
inds <- tibble(id = gl_sec_het_rel_maf@ind.names) %>%
  left_join(meta, by = "id") %>%
  mutate(location = case_when(grepl("oman", pop) ~ "oman",
                              grepl("saudi", pop) ~ "saudi_arabia",
                              grepl("sharjah", pop) ~ "uae",
                              grepl("abudhabi", pop) ~ "uae",
                              grepl("rak", pop) ~ "uae",
                              grepl("dubai", pop) ~ "uae",
                              grepl("srilanka", pop) ~ "srilanka",
                              grepl("bangladesh", pop) ~ "bangladesh")) %>%
  filter(location == "uae") %>%
  sample_n(33-10, replace = F)

uae_remove <- inds$id

gl_sec_het_rel_maf_drop <- gl.drop.ind(gl_sec_het_rel_maf,
                                       ind.list = uae_remove)

gl_sec_het_rel_maf_drop <- gl.filter.monomorphs(gl_sec_het_rel_maf_drop)
gl_sec_het_rel_maf_drop <- gl.recalc.metrics(gl_sec_het_rel_maf_drop) 
nLoc(gl_sec_het_rel_maf_drop)

# edit gl object prior to converting to plink

chr <- tibble(locnames = gl_sec_het_rel_maf_drop@loc.names) %>%
  separate(locnames, c("chr", "pos", "geno"), sep = "-") %>%
  mutate(chr = as.factor(chr),
         pos = as.numeric(pos))

gl_plink <- gl_sec_het_rel_maf_drop

gl_plink@chromosome <- chr$chr

dartR::gl2plink(gl_plink,
                plink_path = "/Users/emilyhumble/software",
                bed_file = T,
                outfile = "gl_sec_het_rhina_drop",
                outpath = "data/out/plink/",
                pos_cM = "0",
                ID_dad = "0",
                ID_mom = "0",
                sex_code = "unknown",
                phen_value = "0",
                verbose = 3)


#~~ Run admixture

system("/Users/emilyhumble/software/admixture_macosx-1.3.0/admixture --cv -B1000 /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/plink/gl_sec_het_rhina_drop.bed 1 -j4 > /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/admixture/gl_sec_het_rhina_drop_admixture_1.out")
system("/Users/emilyhumble/software/admixture_macosx-1.3.0/admixture --cv -B1000 /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/plink/gl_sec_het_rhina_drop.bed 2 -j4 > /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/admixture/gl_sec_het_rhina_drop_admixture_2.out")
system("/Users/emilyhumble/software/admixture_macosx-1.3.0/admixture --cv -B1000 /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/plink/gl_sec_het_rhina_drop.bed 3 -j4 > /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/admixture/gl_sec_het_rhina_drop_admixture_3.out")
system("/Users/emilyhumble/software/admixture_macosx-1.3.0/admixture --cv -B1000 /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/plink/gl_sec_het_rhina_drop.bed 4 -j4 > /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/admixture/gl_sec_het_rhina_drop_admixture_4.out")
system("/Users/emilyhumble/software/admixture_macosx-1.3.0/admixture --cv -B1000 /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/plink/gl_sec_het_rhina_drop.bed 5 -j4 > /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/admixture/gl_sec_het_rhina_drop_admixture_5.out")
system("/Users/emilyhumble/software/admixture_macosx-1.3.0/admixture --cv -B1000 /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/plink/gl_sec_het_rhina_drop.bed 6 -j4 > /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/admixture/gl_sec_het_rhina_drop_admixture_6.out")
system("/Users/emilyhumble/software/admixture_macosx-1.3.0/admixture --cv -B1000 /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/plink/gl_sec_het_rhina_drop.bed 7 -j4 > /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/admixture/gl_sec_het_rhina_drop_admixture_7.out")
system("/Users/emilyhumble/software/admixture_macosx-1.3.0/admixture --cv -B1000 /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/plink/gl_sec_het_rhina_drop.bed 8 -j4 > /Users/emilyhumble/Library/CloudStorage/Dropbox/congen/PROJECTS/IO_elasmo/bowmouth_dart/data/out/admixture/gl_sec_het_rhina_drop_admixture_8.out")


# output files always saved in working directory, move to /admixture

system("mv gl_sec* data/out/admixture/")

system("grep 'CV' data/out/admixture/gl*out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' > data/out/admixture/cv_error_drop")

system("grep 'CV' data/out/admixture/drop/gl*out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' > data/out/admixture/drop/cv_error_drop")


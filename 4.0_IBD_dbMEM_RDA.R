# Script to characterise oversea distance ~ Fst (isolation by distance)
# Using dbMEM by RDA analysis
# Based on mtDNA Fst

library(marmap)
library(data.table)
library(dplyr)
library(adegenet)
library(ggplot2)
library(wesanderson)
library(patchwork)
library(tidyr)
library(hierfstat)
library(scales)
library(geosphere)
library(ggtext)
library(vegan)
library(adespatial)
source("scripts/theme_emily.R")


#~~ Get distance matrix

# Sampling locations

latlong <- fread("data/meta/lat_longs.csv", header = T) %>%
  distinct(pop, .keep_all = T)

# latlong <- latlong %>%
#   filter(pop == "dubai" |
#          pop == "oman_salalah" |
#          pop == "saudiarabia_alqunfudha" |
#          pop == "srilanka" |
#          pop == "bangladesh") %>%
#   mutate(pop = case_when(pop == "oman_salalah" ~ "oman",
#                          pop == "saudiarabia_alqunfudha" ~ "saudiarabia",
#                          pop == "dubai" ~ "uae",
#                          TRUE ~ pop))

#coords <- latlong[,c(6,5)]
coords <- latlong[,c(3,2)]
rownames(coords) <- latlong$pop


# import NOAA bathymetry data
# needs to be high res or else paths go over land

bat <- getNOAA.bathy(lon1 = -180,
                     lon2 = 180,
                     lat1 = -60,
                     lat2 = 30, res = 4, 
                     keep = TRUE, path = "data/marmap/")

# Get depth of points (must all be negative)

get.depth(bat, x = coords$long, y = coords$lat, locator = F)

# Create transition object with a low negative depth constraint to avoid inaccuracies (-10)
# Takes a long time, saved and reloading

# tr <- trans.mat(bat, min.depth = -6, max.depth = NULL)
# save(tr, file = "data/marmap/transition.Rdata")
load("data/marmap/transition.Rdata")

# Compute least cost paths below -10m and plot on map
# Takes a long time, saved and reloading

# lc_path <- lc.dist(tr, coords, res="path")
# save(lc_path, file = "data/marmap/lc_path.Rdata")
load("data/marmap/lc_path.Rdata")

# Distance matrix

# lc_matrix <- lc.dist(tr, coords, res="dist")
# save(lc_matrix, file = "data/marmap/lc_matrix.Rdata")
load("data/marmap/lc_matrix.Rdata")

# check distances make sense and prepare matrix

lc_matrix
geodist <- as.matrix(lc_matrix)
rownames(geodist) <- rownames(coords)
colnames(geodist) <- rownames(geodist)

# Remove AP due to small sample size

#geodist_alf <- geodist_alf[c(1:4,6),c(1:4,6)]

#~~ Prepare genetic info

#~~ COI

pw_fst_coi <- data.frame(OMAN = c(0, 0.30033, 0.07247, 0.52231, 0.07591),
                        SL = c(0.30033, 0, 0.32492, 0.27407, 0.40299),
                        UAE = c(0.07247, 0.32492, 0, 0.49757, 0.04995),
                        BAN = c(0.52231, 0.27407, 0.49757, 0, 0.75000),
                        SA = c(0.07591, 0.40299, 0.04995, 0.75000, 0))

#~~ CR

pw_fst_cr <- data.frame(OMAN = c(0, 0.09568, -0.02693, 0.76252, -0.18736),
                        SL = c(0.09568, 0, 0.13640, 0.06306, 0),
                        UAE = c(-0.02693, 0.13640, 0, 0.65001, -0.15389),
                        BAN = c(0.76252, 0.51020, 0.65001, 0, 0.75000),
                        SA = c(-0.18736, 0, -0.15389, 0.75000, 0))

#~~ Run Mantel test

#~~ Convert to dist object

geodist <- as.dist(geodist)

fstWC_coi <- as.dist(pw_fst_coi)
plot(fstWC_coi ~ geodist)

fstWC_cr <- as.dist(pw_fst_cr)
plot(fstWC_cr ~ geodist)


#~~ Mantel test using all sites
man_coi <- mantel.rtest(geodist, fstWC_coi, nrepet = 1000)
man_coi

man_cr <- mantel.rtest(geodist, fstWC_cr, nrepet = 1000)
man_cr

#~~ dbMEM & RDA analysis

# Transform genetic distances into PCs
fstWC_coi_pc <- prcomp(as.dist(fstWC_coi))
fstWC_cr_pc <- prcomp(as.dist(fstWC_cr))

# Transform geographic distances into dbMEMs
geodist_dbmem <- dbmem(as.dist(geodist))

# RDA
coi_rda <- rda(fstWC_coi_pc$x, geodist_dbmem)
RsquareAdj(coi_rda)
anova(coi_rda, perm=1000)

cr_rda <- rda(fstWC_cr_pc$x, geodist_dbmem)
RsquareAdj(cr_rda)
anova(cr_rda, perm=1000)


#~~ Create labels for plots

fstlab <- expression(italic("F")[ST])

#~~ Create dataframe of distance matrices

df_dist_coi <- data.frame(geodistance=as.vector(geodist),
                          gendistance=as.vector(fstWC_coi))

df_dist_cr <- data.frame(geodistance=as.vector(geodist),
                          gendistance=as.vector(fstWC_cr))

#~~ Plot
col_palette <- c("#F2B705",
                 "#046C9A",
                 "#FD6467",
                 "#9986A5",
                 "#0B775E",
                 "#D67236",
                 "#4F1C2E")


ibd_coi <- ggplot(df_dist_coi) + 
  geom_point(aes(x = geodistance, y = gendistance),
             col = "#046C9A", fill = "#046C9A", # "grey69"
             shape = 21, alpha = 0.5, stroke = 0.1, size = 3) + 
  geom_smooth(aes(x = geodistance, y = gendistance), method = "lm",
              se = TRUE, level= 0.95, alpha = 0.05, size = 1, col = "#046C9A", fill = "#046C9A") +
  xlab("Geographic distance (km)") +
  ylab(fstlab) +
  scale_x_continuous(labels = comma) +
  ggtitle("D") +
  theme_emily() +
  geom_richtext(x = 2500, 
                y = 0.8, 
                label.color = NA,
                label = "Adj *R<sup>2</sup>* = 0.79",
                size = 3,
                color = "#333333")

ibd_coi

ibd_cr <- ggplot(df_dist_cr) + 
  geom_point(aes(x = geodistance, y = gendistance),
             col = "#D67236", fill = "#D67236", # "grey69"
             shape = 21, alpha = 0.5, stroke = 0.1, size = 3) + 
  geom_smooth(aes(x = geodistance, y = gendistance), method = "lm",
              se = TRUE, level= 0.95, alpha = 0.05, size = 1, col = "#D67236", fill = "#D67236") +
  xlab("Geographic distance (km)") +
  ylab(fstlab) +
  scale_x_continuous(labels = comma) +
  ggtitle("B") +
  theme_emily() +
  geom_richtext(x = 2500, 
                y = 1, 
                label.color = NA,
                label = "Adj *R<sup>2</sup>* = 0.55",
                size = 3,
                color = "#333333")

ibd_cr

#saveRDS(ibd_alf, "figs/ibd_alf.RDS")


#~~ Plotting for manuscript

ibd_coi + ibd_cr

#ggsave("figs/IBD.png", ibd_alf + ibd_bir, width = 7, height = 3)

man_coi
r1_coi
p1_coi

#~~ Fst plotting

# Get upper triangle of the correlation matrix

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}


#~~ COI

rownames(pw_fst_coi) <- colnames(pw_fst_coi)
pw_fst_coi <- pw_fst_coi[c(5,1,3,2,4),c(5,1,3,2,4)]
coi_upper_tri <- as.matrix(get_upper_tri(pw_fst_coi))
coi_melted_cormat <- melt(coi_upper_tri, na.rm = TRUE)

coi_fst_plot <- ggplot(coi_melted_cormat, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  #scale_fill_gradient2(name="Fst") +
  scale_fill_gradient(low = "#F8FCFC",
                      high = "#046C9A",
                      name = expression(italic("F")[ST])) + #limit = c(-0.0025,0.160)
  theme_emily() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  coord_fixed() +
  ggtitle("C")
# ggtitle(expression(paste("Reef manta ray ", italic("(Mobula alfredi)"))))

coi_fst_plot

#~~ CR

rownames(pw_fst_cr) <- colnames(pw_fst_cr)
pw_fst_cr <- pw_fst_cr[c(5,1,3,2,4),c(5,1,3,2,4)]
cr_upper_tri <- as.matrix(get_upper_tri(pw_fst_cr))
cr_melted_cormat <- melt(cr_upper_tri, na.rm = TRUE)

cr_fst_plot <- ggplot(cr_melted_cormat, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  #scale_fill_gradient2(name="Fst") +
  scale_fill_gradient(low = "#F8FCFC",
                      high = "#D67236",
                      name = expression(italic("F")[ST])) + #limit = c(-0.0025,0.160)
  theme_emily() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  coord_fixed() +
  ggtitle("A")
# ggtitle(expression(paste("Reef manta ray ", italic("(Mobula alfredi)"))))

cr_fst_plot


#~~ PLOTTING FOR MANUSCRIPT


ggsave("figs/coi_fst_ibd.png", coi_fst_plot + ibd_coi + 
         plot_layout(ncol = 2,
                     widths = c(1, 1),
                     heights = c(1.5, 1)), 
       width = 8, height = 5)

# Calculate MLH across individuals

library(data.table)
library(inbreedR)
library(dplyr)
library(ggplot2)
library(gghalves)
library(patchwork)
source("scripts/theme_emily.R")

# get raw plink from plink files
# calc sMLH


infile <- "data/out/plink/gl_sec_het_rhina"
system(paste0("/Users/emilyhumble/software/plink --bfile ", infile, " --recodeAD --make-bed --allow-extra-chr --debug --out data/out/plink/gl_sec_het_rhina"))


get_MLH_from_plinkraw <- function(file) {
  
  x <- fread(file, colClasses = "character")
  ids <- x$IID
  x <- dplyr::select(x, contains("HET"))
  #row.names(x) <- ids
  NAs <- apply(x, 1, function(x) sum(is.na(x)))
  
  MLH <- as.data.frame(MLH(x))
  MLH$ANIMAL <- ids
  MLH$NAS <- NAs
  colnames(MLH) <- c("MLH", "ANIMAL", "NAs")
  MLH <- MLH
  
}


raw_file <- "data/out/plink/gl_sec_het_rhina.raw"
MLH <- get_MLH_from_plinkraw(raw_file)

meta <- read.csv("data/meta/rhina_metadata.csv")

MLH <- MLH %>%
  left_join(meta, by = c("ANIMAL" = "id")) %>%
  mutate(location = case_when(grepl("oman", pop) ~ "OMAN",
                              grepl("saudi", pop) ~ "SA",
                              grepl("sharjah", pop) ~ "UAE",
                              grepl("abudhabi", pop) ~ "UAE",
                              grepl("rak", pop) ~ "UAE",
                              grepl("dubai", pop) ~ "UAE",
                              grepl("srilanka", pop) ~ "SL",
                              grepl("bangladesh", pop) ~ "BAN")) %>%
  mutate(location = factor(location, levels = c("SA",
                                                "OMAN",
                                                "UAE",
                                                "SL",
                                                "BAN")))

pal <- c("#FD6467",
         "#046C9A",
         "#0B775E",
         "#9986A5",
         "#F2B705",
         "#D67236")


# MLH_filter <- filter(MLH, location != "saudi_arabia" & location != "bangladesh")

het_plot <- ggplot(MLH, aes(as.factor(location), MLH, fill = location)) +
  geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size = 3,
                  transformation_params = list(height = 0, width = 1.3, seed = 1)) +
  geom_half_boxplot(side = "r", outlier.color = NA,
                    width = 0.6, lwd = 0.3, color = "black",
                    alpha = 0.8) +
  theme_emily() +
  scale_fill_manual(values = pal) + 
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  labs(x = "Species", y = "Multi-locus heterozygosity", title = "C")

het_plot

ggsave("figs/het_rhina.png", het_plot, height = 3, width = 4)

#~~ Plotting for manuscript

K2_plot <- readRDS("figs/structure_K2_plot.RDS")
K3_plot <- readRDS("figs/structure_K3_plot.RDS")
K4_plot <- readRDS("figs/structure_K4_plot.RDS")
K5_plot <- readRDS("figs/structure_K5_plot.RDS")
K6_plot <- readRDS("figs/structure_K6_plot.RDS")
K7_plot <- readRDS("figs/structure_K7_plot.RDS")
K8_plot <- readRDS("figs/structure_K8_plot.RDS")

pca_12 <- readRDS("figs/PCA_rhina_12.RDS")


fig_3 <- ((pca_12 / het_plot + plot_layout(guides = 'auto')) | 
            K2_plot / K3_plot / K4_plot / K5_plot / K6_plot) + 
  plot_layout(guides = 'collect', widths = c(1.5, 1))

fig_3

ggsave("figs/fig_3.png", fig_3, height = 7, width = 7)


 #~~ Without Bangladesh

ggplot(filter(MLH, location !="bangladesh"), aes(as.factor(location), MLH, fill = location)) +
  geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size = 3,
                  transformation_params = list(height = 0, width = 1.3, seed = 1)) +
  geom_half_boxplot(side = "r", outlier.color = NA,
                    width = 0.6, lwd = 0.3, color = "black",
                    alpha = 0.8) +
  theme_emily() +
  scale_fill_manual(values = pal) + 
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  labs(x = "Species", y = "Multi-locus heterozygosity")


ID22, ID24, ID58, ID19, ID62, ID54, ID35

ID22 het outlier
ID24 het outlier
ID58 het outlier
ID19
ID62
ID54: 
ID35: high missing values


# Het outliers
# ID25*
# ID33
# ID58*
# ID20
# ID28
# ID29
# ID22*
# ID24*

hist(MLH$MLH)

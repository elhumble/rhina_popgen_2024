# Plotting admixture results

library(purrr)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(wesanderson)
library(patchwork)
library(readxl)
source("scripts/theme_emily.R")


data_path <- "data/out/admixture/drop/"
data_path <- "data/out/admixture/"

files <- dir(data_path, pattern = "*Q$")

admix_data <- tibble(filename = files) %>%
  dplyr::mutate(file_contents = purrr::map(filename,
                                           ~ fread(file.path(data_path, .))))

# add individual names

inds <- fread("data/out/plink/gl_sec_het_rhina.fam") %>%
  dplyr::select(V1, V2) %>%
  rename(ID = V2,
         pop = V1)


# inds <- fread("data/out/plink/gl_sec_het_rhina_drop.fam") %>%
#   dplyr::select(V1, V2) %>%
#   rename(ID = V2,
#          pop = V1)

admix_df <- unnest(admix_data, cols = c(file_contents)) %>%
  mutate(ID = rep(inds$ID, 8)) %>% # number of Ks
  mutate(pop = rep(inds$pop, 8)) %>% # number of Ks
  pivot_longer(cols = c(V1, V2, V3, V4, V5, V6, V7, V8)) %>%
  drop_na() %>%
  mutate(pop = case_when(pop == "UAE" ~ "UAE",
                             pop == "Oman" ~ "OMAN",
                             pop == "Sri_Lanka" ~ "SL",
                             pop == "Saudi_Arabia" ~ "SA",
                         pop == "Bangladesh" ~ "BAN")) %>%
  mutate(pop = fct_relevel(pop, "SA", "OMAN", "UAE", "SL", "BAN"))


K1 <- filter(admix_df, grepl("1.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K2 <- filter(admix_df, grepl("2.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

# Some number for K = 2

K2 %>%
  filter(pop == "OMAN" | pop == "SA") %>%
  filter(name == "V1") %>%
  summarise(mean = mean(value))

K2 %>%
  filter(pop == "SL" | pop == "BAN") %>%
  filter(name == "V2") %>%
  summarise(min = min(value))

K2 %>%
  filter(pop == "UAE") %>%
  filter(name == "V1") %>%
  summarise(mean = mean(value))


K3 <- filter(admix_df, grepl("3.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K4 <- filter(admix_df, grepl("4.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K5 <- filter(admix_df, grepl("5.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K6 <- filter(admix_df, grepl("6.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K7 <- filter(admix_df, grepl("7.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K8 <- filter(admix_df, grepl("8.Q$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))
# Plot

# https://luisdva.github.io/rstats/model-cluster-plots/

col_palette <- c("#F2B705",
                 "#046C9A",
                 "#FD6467",
                 "#9986A5",
                 "#0B775E",
                 "#D67236",
                 "#E2D200",
                 "#35274A",
                 "#4F1C2E")

k1_plot <- ggplot(K1, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~pop, switch = "x", scales = "free", space = "free") +
  #facet_grid(~fct_inorder(pop), switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "B", y = "K=1") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

k2_plot <- ggplot(K2, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~pop, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", y = "K=2") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")


k3_plot <- ggplot(K3, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~pop, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", y = "K=3") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

k4_plot <- ggplot(K4, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~pop, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", y = "K=4") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

k5_plot <- ggplot(K5, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~pop, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", y = "K=5") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

k6_plot <- ggplot(K6, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~pop, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", y = "K=6") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid = element_blank(),
        legend.position = "none")


k7_plot <- ggplot(K7, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~pop, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", y = "K=7") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

k8_plot <- ggplot(K8, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~pop, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", y = "K=8") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        #axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

k1_plot + k2_plot + k3_plot + k4_plot + k5_plot + k6_plot + plot_layout(ncol = 1)


ggsave("figs/admixture_rhina.png",
       k1_plot + k2_plot + k3_plot + k4_plot + k5_plot + k6_plot + plot_layout(ncol = 1),
       width = 6, height = 9)


admixture_plot <- k1_plot + k2_plot + k3_plot + k4_plot + k5_plot + k6_plot + plot_layout(ncol = 1)
saveRDS(admixture_plot, "figs/admixture_plot.RDS")

saveRDS(k1_plot, "figs/admix_K1_plot.RDS")
saveRDS(k2_plot, "figs/admix_K2_plot.RDS")
saveRDS(k3_plot, "figs/admix_K3_plot.RDS")
saveRDS(k4_plot, "figs/admix_K4_plot.RDS")
saveRDS(k5_plot, "figs/admix_K5_plot.RDS")
saveRDS(k6_plot, "figs/admix_K6_plot.RDS")


# Cross validation

cv_err <- fread("data/out/admixture/cv_error")

cv_err <- fread("data/out/admixture/drop/cv_error_drop")

cv_error <- ggplot(cv_err, aes(V1, V2)) +
  geom_point() +
  geom_line() +
  ylab("Cross-validation error") +
  xlab("K") +
  theme_emily() + ggtitle("B")

cv_error


ggsave("figs/admixture_cv.png",
       cv_error,
       width = 5, height = 4)

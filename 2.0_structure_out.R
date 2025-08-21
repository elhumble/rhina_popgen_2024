# Parse and plot STRUCTURE output using mainly pophelper()

library(dplyr)
library(tidyr)
library(stringr)
#install.packages('devtools',dependencies=T)
#library(devtools)
#install_github("royfrancis/pophelper", force = T)
library(pophelper)
library(data.table)
options(scipen=999)
library(wesanderson)
library(forcats)
library(ggplot2)
library(patchwork)
source("scripts/theme_emily.R")


#~~ Metadata and output files

meta <- fread("data/meta/rhina_metadata.csv")

path_to_out <- "data/out/structure"
path_to_out <- "data/out/structure_sub"


sfiles <- list.files(path_to_out, full.names = T, pattern = "_f")
slist <- readQ(files = sfiles, filetype="structure", indlabfromfile = T)
tr1 <- tabulateQ(qlist=slist)
sr1 <- summariseQ(tr1)
em <- evannoMethodStructure(data=sr1)

#~~ Get clumpp input

clumppExport(qlist=slist, exportpath=path_to_out)

# Move to EDDIE to run CLUMPP (can't run CLUMPP on local). See workflow.

# scp -r data/out/structure/pop_K* ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/rhina/
# scp -r data/out/structure_sub/pop_K* ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/rhina/

# On EDDIE qlogin

# for i in {2..8};
# do
# cp /exports/cmvm/eddie/eb/groups/ogden_grp/emily/software/CLUMPP_Linux64.1.1.2/CLUMPP \
# /exports/cmvm/eddie/eb/groups/ogden_grp/emily/rhina/pop_K${i}/.
# done

# for i in {2..8};
# do
# cd pop_K${i}
# ./CLUMPP
# cp pop_K${i}-combined-merged.txt ../.
# rm CLUMPP
# cd ..
# done

# Read in combined CLUMPP output
  
# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/rhina/*-combined-merged.txt data/out/structure/
# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/rhina/*-combined-merged.txt data/out/structure_sub/

# Ks and log likelihoods
  
deltak <- ggplot(em, aes(x = k, y = deltaK)) +
  geom_point(size = 1, col = "grey30") +
  geom_line(size = 1, col = "grey30") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain")) +
  labs(x = "K", y = expression(paste(Delta,italic("K")))) +
  scale_x_continuous(breaks=c(1:6), labels=c(1:6),limits=c(1,6)) +
  theme_emily() +
  labs(title = "A")


loglike <- ggplot(em, aes(x=k, y = elpdmean)) +
  geom_point(size = 1, col = "grey30") + # 1/1.5
  geom_line(size = 1, col = "grey30") +
  geom_errorbar(aes(ymin = elpdmean - elpdsd, ymax= elpdmean + elpdsd), colour="grey30", width=0) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain")) +
  labs(x = "K", y = expression(paste("Ln Pr(",italic("X"),"|",italic("K"),")"))) +
  scale_x_continuous(breaks=c(1:6), labels=c(1:6),limits=c(1,6)) +
  theme_emily() +
  labs(title = "B")


deltak + loglike

ggsave("figs/structure_sub_Ks.png",
       deltak + loglike,
       width = 8, height = 3)

ggsave("figs/structure_Ks.png",
       deltak + loglike,
       width = 8, height = 3)

saveRDS(deltak, "figs/deltak_sub.RDS")
saveRDS(loglike, "figs/loglike_sub.RDS")

#~~ Plotting

data_path <- "data/out/structure"
data_path <- "data/out/structure_sub"

files <- dir(data_path, pattern = ".txt")

data <- tibble(filename = files) %>%
  dplyr::mutate(file_contents = purrr::map(filename,
                                           ~ fread(file.path(data_path, .))))

inds <- fread("data/out/structure/results_job_T1_1_q") %>%
  dplyr::select(V1) %>%
  dplyr::rename(ID = V1)

# Drop last column from each dataframe

data <- data %>%
  mutate(props = purrr::map(data$file_contents, ~ .x[,1:(ncol(.x)-1), drop = F])) %>%
  select(-file_contents)

data_plot <- unnest(data, cols = c(props)) %>%
  mutate(ID = rep(inds$ID, 7)) %>% # number of Ks
  pivot_longer(cols = c(V2, V3, V4, V5, V6, V7, V8, V9)) %>%
  drop_na() %>%
  left_join(meta, by = c("ID" = "id")) %>%
  mutate(pop = case_when(grepl("oman", pop) ~ "OMAN",
                              grepl("saudi", pop) ~ "SA",
                              grepl("sharjah", pop) ~ "UAE",
                              grepl("abudhabi", pop) ~ "UAE",
                              grepl("rak", pop) ~ "UAE",
                              grepl("dubai", pop) ~ "UAE",
                              grepl("srilanka", pop) ~ "SL",
                              grepl("bangladesh", pop) ~ "BAN"))

K2 <- filter(data_plot, grepl("K2", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

# Some numbers for K = 2

K2 %>%
  filter(pop == "OMAN" | pop == "SA") %>%
  filter(name == "V3") %>%
  summarise(mean = mean(value))

K2 %>%
  filter(pop == "SL" | pop == "BAN") %>%
  filter(name == "V2") %>%
  summarise(min = min(value))

K2 %>%
  filter(pop == "UAE") %>%
  filter(name == "V2") %>%
  summarise(mean = mean(value))

K2 %>%
  filter(pop == "UAE") %>%
  filter(name == "V3") %>%
  summarise(mean = mean(value))

#~~~

K3 <- filter(data_plot, grepl("K3", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K4 <- filter(data_plot, grepl("K4", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K5 <- filter(data_plot, grepl("K5", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K6 <- filter(data_plot, grepl("K6", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K7 <- filter(data_plot, grepl("K7", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K8 <- filter(data_plot, grepl("K8", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

col_palette <- c("#F2B705",
                 "#046C9A",
                 "#FD6467",
                 "#9986A5",
                 "#0B775E",
                 "#D67236",
                 "#E2D200",
                 "#35274A",
                 "#4F1C2E")


k2_plot <- ggplot(K2, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(pop), switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "B", y = "K=2") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(size=10),
        panel.grid = element_blank(),
        legend.position = "none")

k3_plot <- ggplot(K3, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(pop), switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", y = "K=3") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(size=10),
        panel.grid = element_blank(),
        legend.position = "none")

k4_plot <- ggplot(K4, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(pop), switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", y = "K=4") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size=10),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

k5_plot <- ggplot(K5, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(pop), switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", y = "K=5") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size=10),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")


k6_plot <- ggplot(K6, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(pop), switch = "x", scales = "free", space = "free") +
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
  facet_grid(~fct_inorder(pop), switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", y = "K=7") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size=10),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

k8_plot <- ggplot(K8, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(pop), switch = "x", scales = "free", space = "free",
             labeller = label_wrap_gen(width = 10, multi_line = T)) +
  theme_minimal() + labs(x = "Individuals", y = "K=8") +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size=10),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")


k2_plot + k3_plot + k4_plot + k5_plot + k6_plot + k7_plot + k8_plot + plot_layout(ncol = 1)

# Updated figure for manuscript

k2_panel <- ggplot(K2, aes(value, factor(ID), fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(rows = vars(fct_inorder(pop)), switch = "y", scales = "free", space = "free") +
  theme_minimal() + labs(y = "Individuals", title = "B", x = "K=2") +
  scale_fill_manual(values = col_palette) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(size=10),
        panel.grid = element_blank(),
        legend.position = "none")

# Save figure files

saveRDS(k2_plot, "figs/structure_K2_plot.RDS")
saveRDS(k3_plot, "figs/structure_K3_plot.RDS")
saveRDS(k4_plot, "figs/structure_K4_plot.RDS")
saveRDS(k5_plot, "figs/structure_K5_plot.RDS")
saveRDS(k6_plot, "figs/structure_K6_plot.RDS")
saveRDS(k7_plot, "figs/structure_K7_plot.RDS")
saveRDS(k8_plot, "figs/structure_K8_plot.RDS")
saveRDS(k2_panel, "figs/structure_K2_panel.RDS")


saveRDS(k2_plot, "figs/structure_sub_K2_plot.RDS")
saveRDS(k3_plot, "figs/structure_sub_K3_plot.RDS")
saveRDS(k4_plot, "figs/structure_sub_K4_plot.RDS")
saveRDS(k5_plot, "figs/structure_sub_K5_plot.RDS")
saveRDS(k6_plot, "figs/structure_sub_K6_plot.RDS")
saveRDS(k7_plot, "figs/structure_sub_K7_plot.RDS")
saveRDS(k8_plot, "figs/structure_sub_K8_plot.RDS")


ggsave("figs/structure_rhina.png",
       k2_plot + k3_plot + k4_plot + k5_plot + k6_plot + k7_plot + k8_plot + plot_layout(ncol = 1),
       width = 6, height = 10)


ggsave("figs/structure_sub_rhina.png",
       k2_plot + k3_plot + k4_plot + k5_plot + k6_plot + k7_plot + k8_plot + plot_layout(ncol = 1),
       width = 6, height = 10)


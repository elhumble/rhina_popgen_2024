# Script to plot sample map

library(ggplot2)
library(dplyr)
library(sf)
library(tidyterra)
library(terra)
library(raster)
library(patchwork)
library(wesanderson)
source("scripts/theme_emily.R")

#~~ Read in map data

world.map <- vect("data/map/TM_WORLD_BORDERS_SIMPL-0.3.shp")

# Crop to avoid weird lines
#world.map <- terra::crop(world.map, extent(-179.9, 179.9, -90, 84))

# Make sf object
world.map <- st_as_sf(world.map)

# Remove Antarctica
world.map <- world.map[world.map$NAME != "Antarctica",]

#~~ Read in distribution data from the IUCN

rhina <- vect("data/map/redlist_species_data_fe8e4a1f-c967-4846-a54d-592bd3317bfb/data_0.shp")

rhina_extant <- rhina[rhina$LEGEND=='Extant (resident)',]
rhina_pos_extant <- rhina[rhina$LEGEND=='Possibly Extant (resident)',]

rhina_extant <- st_as_sf(rhina_extant)
rhina_pos_extant <- st_as_sf(rhina_pos_extant)

#~~ Read in sampling location files

sample_sites <- read.csv("data/meta/rhina_metadata_lat_long.csv") %>%
  mutate(Country = case_when(pop == "abudhabi" | pop == "dubai" | pop == "rak" |
                                pop == "sharjah" ~ "UAE",
                              pop == "oman_salalah" | pop == "oman_shinas" | pop == "oman_sohar" ~ "Oman",
                              pop == "saudiarabia_alqunfudha" | pop == "saudiarabia_jeddah" | 
                                pop == "saudiarabia_jizan" ~ "Saudi Arabia",
                              pop == "bangladesh" ~ "Bangladesh",
                              pop == "srilanka" ~ "Sri Lanka",
                             pop == "taiwan" ~ "Taiwan"))


# Summarise number of samples per site

rhina_samples <- sample_sites %>%
  dplyr::group_by(pop, Country, lat, long) %>%
  dplyr::summarise(n = n()) %>%
  mutate(lon_samples = long,
         lat_samples = lat)

sample_sites %>%
  dplyr::group_by(Country) %>%
  dplyr::summarise(n = n())

#~~ Plot map

pal <- c("#F2B705",
         "#046C9A",
         "#0B775E",
         "#9986A5",
         "#D67236",
         "#FD6467",
         "#4F1C2E")

pal <- c("#F2B705",
         "#046C9A",
         "#FD6467",
         "#9986A5",
         "#D67236",
         "#0B775E",
         "#4F1C2E")

m_rhina <- ggplot() + 
  geom_spatvector(data = world.map, colour = "black", fill = "#F8F1EB", size = 0.1) +
 # geom_spatvector(data = rhina_extant, col = NA, fill = "#7E8AA2", alpha = 1) +
  geom_spatvector(data = rhina_extant, col = NA, fill = "#526180", alpha = 1) +
  geom_point(aes(x = lon_samples, y = lat_samples, color = Country, size = n), 
             data = rhina_samples, alpha = 0.9) + # 3, 2
  scale_colour_manual(values = pal) + 
  xlim(c(-15, 160)) + 
  ylim(c(-45, 50)) +
  # theme_emily() +
  #guides(color = "none") +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#e4eff7"),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.04),
        plot.subtitle = element_text(hjust = 0.08)) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  scale_size(range = c(3,9))
  #ggtitle("C")
  #ggtitle(label = "C",
   #       subtitle = expression(paste("Bowmouth Guitarfish ", italic("(Rhina ancylostoma)"))))

ggsave("figs/map_rhina.png", m_rhina, width = 10, height = 5)
ggsave("figs/map_rhina.pdf", m_rhina, width = 10, height = 5)

 # Distribution only

m_dist <- ggplot() + 
  geom_spatvector(data = world.map, colour = "black", fill = "#F8F1EB", size = 0.1) +
  geom_spatvector(data = rhina_extant, col = NA, fill = "#7E8AA2", alpha = 0.9) +
  xlim(c(-15, 160)) + 
  ylim(c(-45, 50)) +
  # theme_emily() +
  #guides(color = "none") +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#e4eff7"),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.04),
        plot.subtitle = element_text(hjust = 0.08)) +
  scale_size(range = c(3,9))
#ggtitle("C")
#ggtitle(label = "C",
#       subtitle = expr

ggsave("figs/map_dist.png", m_dist, width = 11, height = 6)

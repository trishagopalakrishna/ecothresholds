library(tidyverse)
library(here)
library(ggplot2)
library(ggpubr)
library(sf)
library(terra)
library(tmap)
library(tmaptools)
library(RColorBrewer)
library(viridis)

input_filepath <- here("Outputs", "TrendsResults", "Figures")

byveg<- read_rds(here(input_filepath, "agg_trajectories_vegformations_boxplot.Rdata"))
byclimatezone <- read_rds(here(input_filepath, "agg_trajectories_climatezones.Rdata"))
byfloristicregions <- read_rds(here(input_filepath, "agg_trajectories_floristicregions.Rdata"))
byfireregions <- read_rds(here(input_filepath, "agg_trajectories_fireregions.Rdata"))

#1. Vegetation formation
vegformation_palette<- c("#638b66", "#b66353", "#fbb04e" )

y_boxplot <- ggplot(byveg, aes(y = value, x = Vegformation)) +
  geom_boxplot(aes(fill = Vegformation)) + 
  facet_wrap(.~TrajType) +
  theme_classic(base_size = 14) +
  theme(legend.title=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+ 
  scale_fill_manual(values = vegformation_palette) +
  xlab("Vegetation Formation") + ylab("Proportion of 1km pixel")

aggregate_trajectory_palette <- c("#663300", "#61A36A", "#7C7083")

#2. Climate zones
climate_zones <- ggplot(byclimatezone, aes(x = factor(Region), y= area_sqkm/1000, fill = factor(TrajectoryType))) + 
  geom_bar(stat = "identity") + theme_classic() +
  scale_fill_manual(values =aggregate_trajectory_palette)+
  theme_classic(base_size = 14) +
  theme(legend.position="none")+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) + ylab ("Area (*1000 sqkm)") 

#3. Floristic zones
floristic_zones <- ggplot(byfloristicregions, aes(x = factor(Region, 
                          levels = c("Central-East", "Central-West", "North-East", "North-West", 
                                     "South", "South-East", "South-West", "external North")), 
                          y= area_sqkm/1000, fill =factor(TrajectoryType))) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values =aggregate_trajectory_palette)+
  theme_classic(base_size = 14) +
  theme(legend.position="none")+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust=0.3)) + ylab ("Area (*1000 sqkm)") 

#4. Fire regions
fire_zones <- ggplot(byfireregions, aes(x = reorder(factor(Region, 
                                        levels = c("LS_P4", "HFA_P4", "LFA_P4")), -(area_sqkm/1000)), 
                             y= area_sqkm/1000, fill =factor(TrajectoryType))) + 
  geom_bar(stat = "identity")  +
  scale_fill_manual(values =aggregate_trajectory_palette)+
  theme_classic(base_size = 14) +
  theme(legend.position="none")+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) + ylab ("Area (*1000 sqkm)") 




fig2 <- ggpubr::ggarrange(plotlist = list(y_boxplot, climate_zones, floristic_zones, fire_zones),
                          ncol = 2, nrow = 2, labels = c("a)", "b)", "c)", "d)"))
ggsave(here(input_filepath, "fig2.png"), fig2 ,dpi = 700, height = 25, width = 25, units = "cm")

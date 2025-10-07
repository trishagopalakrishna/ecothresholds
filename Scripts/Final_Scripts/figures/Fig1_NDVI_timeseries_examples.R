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

############################## Fig 1a - Map of NDVI based greening, browning and no trend trajectories + examples of types of trajectories
#1. Map 
d_trans<- st_read(here ("Data", "Cattelanetal_Clustering", "cerrado_climate_zones.shp"))
d_trans<- st_transform(d_trans, crs = 4326)
d_trans<- d_trans %>% mutate(area_km= terra::expanse(vect(d_trans), unit="km"))
cerrado<- d_trans %>% st_union()

monthly_ts_swin11_results <- rast(here("Outputs", "TrendsResults", "aggregate_trajectory_results", "monthly", "reclass_monthly_ndvi_swin11.tif"))

aggregate_trajectory_palette <- c("#663300", "#61A36A", "#7C7083")

aggregate_mapping_function <- function (trendresults_raster){
  x_map<- tm_shape(cerrado)+ tm_borders()+
    tm_shape (trendresults_raster) + 
    tm_raster(col.scale = tm_scale_categorical(n=3, values = aggregate_trajectory_palette), 
              col.legend = tm_legend_hide())
  x_map
}

fig1_map <- aggregate_mapping_function(monthly_ts_swin11_results)
output_file_path <- here("Outputs", "TrendsResults", "aggregate_trajectory_results")
tmap_save(fig1_map, paste0(output_file_path, "/", "map1km_ndvi_aggtrajectories_monthly_swin11.png"),
          height = 80, width = 80, units = "cm", dpi=700)

#2. Examples of time series of different types of trajectories
ts_data <- read.csv(here("Outputs", "TrajectoryPlotting", "ndvi_evi_ts_trajectories.csv"))

lp_fig<- function (df, title) {
  df %>%  
  ggplot(aes(x = as.factor(time), y= NDVI_trend, group=1))+
  geom_line(linetype="solid")+
  theme_classic()+
    scale_x_discrete(breaks = c("2002_01", "2008_01", "2014_01", "2021_01"), labels = c("2002", "2008", "2014", "2021"))+
  theme (axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.title.x = element_blank(),
         axis.text.x = element_text(angle = 45, hjust=1),
         plot.title = element_text(hjust = 0.5))+
  ylab("NDVI decomposed monthly value") + ggtitle(title)+
   scale_alpha(guide = 'none')
}

lin_inc <- ts_data %>% filter(NDVI_ID==200529) %>% dplyr::select(c(time, NDVI_trend)) %>% mutate(Time2=time) %>%
  separate_wider_delim(Time2, "_", names=c("Year", "Month"))
fig_lin_inc <- lp_fig(lin_inc, "Linear increase")
ggsave(here("Outputs", "TrajectoryPlotting", "maintext_NDVI_lineplots", "lin_inc.png"),
       fig_lin_inc, dpi = 700, height = 7, width = 7, units = "cm")

lin_dec <- ts_data %>% filter(NDVI_ID==3612136) %>% dplyr::select(c(time, NDVI_trend)) %>% mutate(Time2=time) %>%
  separate_wider_delim(Time2, "_", names=c("Year", "Month"))
fig_lin_dec <- lp_fig(lin_dec, "Linear decrease")
ggsave(here("Outputs", "TrajectoryPlotting", "maintext_NDVI_lineplots", "lin_dec.png"),
       fig_lin_dec, dpi = 700, height = 7, width = 7, units = "cm")

step_inc <- ts_data %>% filter(NDVI_ID==1293359) %>% dplyr::select(c(time, NDVI_trend)) %>% mutate(Time2=time) %>%
  separate_wider_delim(Time2, "_", names=c("Year", "Month"))
fig_step_inc <- lp_fig(step_inc, "Abrupt positive")
ggsave(here("Outputs", "TrajectoryPlotting", "maintext_NDVI_lineplots", "step_inc.png"),
       fig_step_inc, dpi = 700, height = 7, width = 7, units = "cm")

step_dec <- ts_data %>% filter(NDVI_ID==1280327) %>% dplyr::select(c(time, NDVI_trend)) %>% mutate(Time2=time) %>%
  separate_wider_delim(Time2, "_", names=c("Year", "Month"))
fig_step_dec <- lp_fig(step_dec, "Abrupt negative")
ggsave(here("Outputs", "TrajectoryPlotting", "maintext_NDVI_lineplots", "step_dec.png"),
       fig_step_dec, dpi = 700, height = 7, width = 7, units = "cm")

no_trend <- ts_data %>% filter(NDVI_ID==8827783) %>% dplyr::select(c(time, NDVI_trend)) %>% mutate(Time2=time) %>%
  separate_wider_delim(Time2, "_", names=c("Year", "Month"))
fig_no_trend <- lp_fig(no_trend, "No trend")
ggsave(here("Outputs", "TrajectoryPlotting", "maintext_NDVI_lineplots", "no_trend.png"),
       fig_no_trend, dpi = 700, height = 7, width = 7, units = "cm")


############################## Fig 1b - Characterization of browning, greening and no trends in terms of vegetation formation
input_filepath <- here("Outputs", "TrendsResults", "Figures")

byveg<- read_rds(here(input_filepath, "agg_trajectories_vegformations_boxplot.Rdata"))

#1. Vegetation formation
vegformation_palette<- c("#638b66", "#b66353", "#fbb04e" )

y_boxplot <- ggplot(byveg, aes(y = value, x = Vegformation)) +
  geom_boxplot(aes(fill = Vegformation), outlier.shape = NA) + 
  facet_wrap(.~TrajType) +
  theme_classic(base_size = 14) +
  theme(legend.title=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+ 
  scale_fill_manual(values = vegformation_palette) +
  xlab("Vegetation Formation") + ylab("Proportion of 1km pixel")
ggsave(here("Outputs","TrendsResults", "Figures", "vegcomposition.png"), y_boxplot,
       dpi = 700, height = 15, width = 15, units = "cm")


############################## Fig 1c - Composition of greening, browning, no trend in terms of trajectory types
trajectorycomposition <- read_rds(here(input_filepath, "trajectorytypes_monthly_swin11.rdata"))

Set1palette_barplot<- c("#994F00",   "#88CCEE" ,  "#332288", "#117733", "lightgrey", "#AA4499", "#882255", "#44AA99", "#CC6677") 

monthly_ndvi_seperate_plot <- trajectorycomposition %>%
  ggplot(aes(x = factor(agg_traj), y= count/1000, fill = factor(seperate_traj))) + 
  geom_bar(stat = "identity")  + 
  theme_classic(base_size = 14) +
  scale_fill_manual(values =Set1palette_barplot)+
  theme(legend.title=element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust=0.7),
        axis.title.x = element_blank(),
        legend.title=element_blank(),
        axis.ticks.x = element_blank()) + ylab ("Area (*1000 sqkm)")
ggsave(here("Outputs","TrendsResults", "Figures", "disaggregate_traj.png"), monthly_ndvi_seperate_plot,
       dpi = 700, height = 15, width = 15, units = "cm")


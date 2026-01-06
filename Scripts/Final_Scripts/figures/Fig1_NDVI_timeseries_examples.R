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
              col.legend = tm_legend(show = F)) +
    tm_layout(frame = FALSE)
  x_map
}

fig1_map <- aggregate_mapping_function(monthly_ts_swin11_results)
output_file_path <- here("Outputs", "TrendsResults", "aggregate_trajectory_results")
tmap_save(fig1_map, paste0(output_file_path, "/", "map1km_ndvi_aggtrajectories_monthly_swin11_2.png"),
          height = 12, width = 12, units = "cm", dpi= 400)

#2. Examples of time series of different types of trajectories
ts_data <- read.csv(here("Outputs", "TrajectoryPlotting", "ndvi_evi_ts_trajectories.csv"))

lp_fig<- function (df, title) {
  df %>%  
  ggplot(aes(x = as.factor(time), y= NDVI_trend, group=1))+
  geom_line(linetype="solid")+
  theme_classic(base_size = 10)+
    scale_x_discrete(breaks = c("2002_01", "2008_01", "2014_01", "2021_01"), labels = c("2002", "2008", "2014", "2021"))+
  theme (axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.title.x = element_blank(),
         axis.text.x = element_text(angle = 45, hjust=1, colour = "black"),
         plot.title = element_text(size = 8,hjust = 0.5))+
  ylab("NDVI deseasoned\nmonthly value") + ggtitle(title)+
   scale_alpha(guide = 'none')
}

lin_inc <- ts_data %>% filter(NDVI_ID==200529) %>% dplyr::select(c(time, NDVI_trend)) %>% mutate(Time2=time) %>%
  separate_wider_delim(Time2, "_", names=c("Year", "Month"))
fig_lin_inc <- lp_fig(lin_inc, "Linear greening")
ggsave(here("Outputs", "TrajectoryPlotting", "maintext_NDVI_lineplots", "lin_inc.png"),
       fig_lin_inc, dpi = 300, height = 3, width = 3, units = "cm")

lin_dec <- ts_data %>% filter(NDVI_ID==3612136) %>% dplyr::select(c(time, NDVI_trend)) %>% mutate(Time2=time) %>%
  separate_wider_delim(Time2, "_", names=c("Year", "Month"))
fig_lin_dec <- lp_fig(lin_dec, "Linear browning")
ggsave(here("Outputs", "TrajectoryPlotting", "maintext_NDVI_lineplots", "lin_dec.png"),
       fig_lin_dec, dpi = 300, height = 3, width = 3, units = "cm")

step_inc <- ts_data %>% filter(NDVI_ID==1293359) %>% dplyr::select(c(time, NDVI_trend)) %>% mutate(Time2=time) %>%
  separate_wider_delim(Time2, "_", names=c("Year", "Month"))
fig_step_inc <- lp_fig(step_inc, "Rapid greening")
ggsave(here("Outputs", "TrajectoryPlotting", "maintext_NDVI_lineplots", "step_inc.png"),
       fig_step_inc, dpi = 300, height = 3, width = 3, units = "cm")

step_dec <- ts_data %>% filter(NDVI_ID==1280327) %>% dplyr::select(c(time, NDVI_trend)) %>% mutate(Time2=time) %>%
  separate_wider_delim(Time2, "_", names=c("Year", "Month"))
fig_step_dec <- lp_fig(step_dec, "Rapid browning")
ggsave(here("Outputs", "TrajectoryPlotting", "maintext_NDVI_lineplots", "step_dec.png"),
       fig_step_dec, dpi = 300, height = 3, width = 3, units = "cm")

no_trend <- ts_data %>% filter(NDVI_ID==8827783) %>% dplyr::select(c(time, NDVI_trend)) %>% mutate(Time2=time) %>%
  separate_wider_delim(Time2, "_", names=c("Year", "Month"))
fig_no_trend <- lp_fig(no_trend, "No transition")
ggsave(here("Outputs", "TrajectoryPlotting", "maintext_NDVI_lineplots", "no_trend.png"),
       fig_no_trend, dpi = 300, height = 3, width = 3, units = "cm")


############################## Fig 1b - Characterization of browning, greening and no trends in terms of vegetation formation
input_filepath <- here("Outputs", "TrendsResults", "Figures")

byveg<- read_rds(here(input_filepath, "agg_trajectories_vegformations_boxplot.Rdata"))

#1. Vegetation formation
#vegformation_palette<- c("#638b66", "#b66353", "#fbb04e" )

y_boxplot <- ggplot(byveg, aes(y = value, x = TrajType)) +
  geom_boxplot(aes(fill = TrajType), outlier.shape = NA) + 
  facet_wrap(.~Vegformation) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+ 
  scale_fill_manual(values = aggregate_trajectory_palette) +
  xlab("Transitions") + ylab("Proportion of 1km pixel\ncovered by each vegetation formation")


############################## Fig 1c - Composition of greening, browning, no trend in terms of trajectory types
trajectorycomposition <- read_rds(here(input_filepath, "trajectorytypes_monthly_swin11.rdata"))
trajectorycomposition[1,1] <- "Rapid browning"
trajectorycomposition[2,1] <- "Rapid greening"
trajectorycomposition[5,1] <- "No transitions"
trajectorycomposition[6,1] <- "Quadratic browning (+)"
trajectorycomposition[7,1] <- "Quadratic browning (-)"
trajectorycomposition[8,1] <- "Quadratic greening (+)"
trajectorycomposition[9,1] <- "Quadratic greening (-)"

#Set1palette_barplot<- c("#994F00","#117733", "#bf812d", "#41AB5D", "lightgrey", "#dfc27d", "#f6e8c3", "#74C476","#C7E9C0")
Set1palette_barplot<- c("#663300","#61A36A", "#bf812d", "#41AB5D" , "#7C7083", "#dfc27d", "#f6e8c3", "#74C476","#C7E9C0")

trajectorycomposition <- trajectorycomposition %>%
  mutate(seperate_traj = factor(seperate_traj, 
                                   levels = c("Rapid browning", "Rapid greening",
                                          "Linear browning", "Linear greening",
                                          "No transitions", "Quadratic browning (+)",
                                          "Quadratic browning (-)", "Quadratic greening (+)",
                                          "Quadratic greening (-)")))

monthly_ndvi_seperate_plot <- trajectorycomposition %>%
  ggplot(aes(x = factor(agg_traj), y= count/1000, fill = factor(seperate_traj))) + 
  geom_bar(stat = "identity")  + 
  theme_classic(base_size = 12) +
  scale_fill_manual(values =Set1palette_barplot) +
  theme(legend.title=element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust=0.4),
        axis.title.x = element_blank(),
        legend.title=element_blank(),
        axis.ticks.x = element_blank()) + ylab ("Proportional area\n (*1000 sqkm)")


x <- ggarrange(monthly_ndvi_seperate_plot, y_boxplot, ncol = 2, nrow = 1, labels = c("b)", "c)"))
ggsave(here("Outputs","TrendsResults", "Figures", "transition_veg_composition2.png"), x,
       dpi = 700, height = 8, width = 18, units = "cm")

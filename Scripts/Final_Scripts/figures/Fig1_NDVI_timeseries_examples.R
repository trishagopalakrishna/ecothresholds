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

#1. Time series input data of NDVI (monthly time step swin=11)
ts_data <- read.csv(here("Outputs", "TrajectoryPlotting", "ndvi_evi_ts_trajectories.csv"))

#2. Line plot function
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


#3.Line plots
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

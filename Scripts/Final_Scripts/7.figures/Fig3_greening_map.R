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
library(mgcv)
library(boot)
library(scales)

#1.Input needed data and models
analyses_df_noNA_nooutliers <- read_rds(here("Scripts", "Final_Scripts", "6.analyses", "analyses_df_noNA_nooutliers.rdata"))

greening_gam <- read_rds(here("Scripts", "Final_Scripts", "6.analyses", "onemodel_inc.rds"))

#2. Greening raster created from data
prediction_raster <- function(gam_model){
  pred <- mgcv::predict.bam(gam_model, newdata = analyses_df_noNA_nooutliers, type="terms") 
  pred <- as.data.frame(pred)
  
  pred_without_xy <- pred %>% dplyr::select(-c("s(x,y)"))
  
  x_effect <- rowSums(pred_without_xy)
  x_effect_inv <- inv.logit(x_effect + gam_model$coefficients[1])
  x_effect_inv <- as.data.frame(x_effect_inv)
  x_effect_inv <- x_effect_inv %>% bind_cols(analyses_df_noNA_nooutliers %>% 
                                         dplyr::select(c(x, y)))
  names(x_effect_inv)[1]<- "main_effect"
  x_effect_inv <- x_effect_inv %>% dplyr::select(x,y, main_effect)
  x_rast <- rast(as.data.frame(x_effect_inv %>% drop_na()), type="xyz")
  
  x_rast
}

Sys.time(); greening_rast <- prediction_raster(greening_gam); Sys.time()
writeRaster(greening_rast, here ("Outputs","AnalysesResults", "GAMResults", "FullModel_Rasters", "full_model_greening.tif"))

#2  Maps of greening and browning 
greening_rast <- rast(here ("Outputs","AnalysesResults", "GAMResults", "FullModel_Rasters", "full_model_greening.tif"))

d_trans<- st_read(here ("Data", "Cattelanetal_Clustering", "cerrado_climate_zones.shp"))
d_trans<- st_transform(d_trans, crs = 4326)
d_trans<- d_trans %>% mutate(area_km= terra::expanse(vect(d_trans), unit="km"))
cerrado<- d_trans %>% st_union()

brazil_states <- st_read(here("Data", "Admin", "BR_UF_2020", "BR_UF_2020.shp"))
brazil_states <- st_transform(brazil_states, crs= st_crs(cerrado))
cerrado_states <- st_intersection(brazil_states, cerrado)
cerrado_states <- cerrado_states %>% filter(NM_UF!="Rondônia" & NM_UF!= "Pará")

aggregate_mapping_function <- function (partialeffect_raster, greens_or_brown_values){
  x_map<- tm_shape(cerrado)+ tm_borders()+
    tm_polygons(fill.scale = tm_scale_categorical(n=1, values = "grey"))+
    tm_shape (partialeffect_raster) + 
    tm_raster(col.scale = tm_scale_continuous(values = greens_or_brown_values),
              col.legend = tm_legend(show= TRUE, title = "Probability", reverse = T, frame.color = NA)) + 
    tm_shape(cerrado_states) + tm_borders(col = "black", lwd = 1.5)+
    tm_text("SIGLA_UF", size = 1, fontface = "bold", col = "#DC267F" )+
    tm_layout(frame = FALSE)
  x_map
}

Sys.time(); map_greening <- aggregate_mapping_function(greening_rast, "greens"); Sys.time()

tmap_save(map_greening, here("Outputs", "TrendsResults", "Figures", "Fig4", "map_green_only.png"))

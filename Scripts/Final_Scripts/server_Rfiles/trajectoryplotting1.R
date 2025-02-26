##Extraction of time series for each index for different trajectories
library(tidyverse)
library(terra)
library(here)
library(sf)
library(ggpubr)
library(tmap)
library(RColorBrewer)


#Data needed in function
d_trans<- st_read(here ("Data", "Cattelanetal_Clustering", "cerrado_climate_zones.shp"))
d_trans<- st_transform(d_trans, crs = 4326)
d_trans<- d_trans %>% mutate(area_km= terra::expanse(vect(d_trans), unit="km"))
cerrado<- d_trans %>% st_union()

evi<- rast(here("Data", "Indices", "AnisoEVI_Cerrado_GEE", "mosaiced_anisoEVI.tif"))
Sys.time(); evi<- terra::crop(evi[[262]], vect(cerrado)); Sys.time()
Sys.time(); evi<- terra::mask(evi, vect(cerrado)); Sys.time()

threshold_mapbiomas_mask_20<-rast(here("Data", "Masks_to_define_pixelsofinterest", "mapbiomas_anthropic_mask", "thresholded_mapbiomas_mask_20.tif"))

#Df prep for each climate zone seperately, making raster and masking
step1<- function (final_results_df, index_name, climate_zone){
  #df prep
  final_results_df<- final_results_df %>% 
    dplyr::select(c(cell,x,y, model_order, shape_class, trend_name, loc_brk, climate_zone))
  
  final_results_df<- final_results_df %>% ungroup() %>%
    dplyr::select(c(cell, x, y, model_order, shape_class, trend_name, climate_zone)) %>%
    mutate(pixelvalue= case_when(model_order=="Lin" & shape_class== "decrease_constant"& is.na(trend_name)~1,
                                 model_order=="Lin" & shape_class== "increase_constant" & is.na(trend_name)~2,
                                 model_order=="Lin" & shape_class== "stable_constant" & is.na(trend_name)~3,
                                 model_order=="Null" & shape_class== "stable_constant" & is.na(trend_name)~3,
                                 model_order=="Quad" & shape_class== "stable_concave" & is.na(trend_name)~3,
                                 model_order=="Quad" & shape_class== "stable_convex" & is.na(trend_name)~3,
                                 model_order=="Step" & is.na(shape_class) & trend_name == "decrease" ~4,
                                 model_order=="Step" & is.na(shape_class) & trend_name == "increase" ~5,
                                 model_order=="Quad" & shape_class=="decrease_accelerated" &is.na(trend_name)~6,
                                 model_order=="Quad" & shape_class=="decrease_decelerated" &is.na(trend_name)~7)) %>%
    dplyr::select(-c(model_order, shape_class,trend_name)) #note Linear & Null stable is same category
  
  #map making
  unique_traj<- unique(final_results_df$pixelvalue)
  
  for (i in 1:length(unique_traj)){
    selected_traj<- final_results_df %>% filter(pixelvalue == i)
    selected_traj_vector<- terra::vect(selected_traj, geom=c("x", "y"), crs="epsg:4326")
    Sys.time(); selected_traj_raster<- terra::rasterize(selected_traj_vector, evi, "cell", fun="max"); Sys.time()
    masked_selected_traj_raster<- terra::mask(selected_traj_raster, threshold_mapbiomas_mask_20)
    writeRaster(masked_selected_traj_raster, here("Outputs", "TrajectoryPlotting", index_name, climate_zone, paste0(i, "_traj_cellnumber.tif")))
  }; Sys.time()
  
  #Matching to df
  cellraster_files<- list.files(path = here("Outputs", "TrajectoryPlotting", index_name, climate_zone), pattern =".tif", all.files= TRUE, full.names = TRUE)
  cellraster_files<- lapply (cellraster_files, rast)  
  
  pixels_studyarea_eachtraj <- list()
  for ( i in 1:length(cellraster_files)){
    raster_to_df <- terra:: as.data.frame(cellraster_files[[i]], xy=TRUE)
    names(raster_to_df)[3]<- "cell"
    pixels_studyarea_eachtraj[[i]] <- left_join (raster_to_df,  final_results_df, by= c("cell"= "cell"))
  }
  
  #Sampling 10 random pixels in each trajctory ie each of the dfs in the list of dfs 
  #from previous step
  randomsample_studyarea_eachtraj <- list()
  for (i in 1:length(pixels_studyarea_eachtraj)){
    randomsample_studyarea_eachtraj[[i]]<- pixels_studyarea_eachtraj[[i]] %>% 
                                            ungroup() %>%
                                              slice_sample(n = 10)
  }
  
  randomsample_studyarea_eachtraj<- bind_rows(randomsample_studyarea_eachtraj)
  randomsample_studyarea_eachtraj<- randomsample_studyarea_eachtraj %>% 
      select(x.x, y.x, cell, climate_zone, pixelvalue) %>% rename("x"="x.x", "y"="y.x")
  write.csv(randomsample_studyarea_eachtraj, here("Outputs", "TrajectoryPlotting", index_name, climate_zone, "sampled_df.csv" ))
}

#Application of above function
#anisoEVI
am_trend_southern<- read_rds(here("Outputs", "Final_TrajShape_Models", "Annual", "finalshape_annuals11_southern.rds"))
am_trend_southern<- am_trend_southern %>% mutate(climate_zone = "southern")
am_trend_eastern<- read_rds(here("Outputs", "Final_TrajShape_Models", "Annual", "finalshape_annuals11_eastern.rds"))
am_trend_eastern<- am_trend_eastern %>% mutate(climate_zone = "eastern")
am_trend_central<- read_rds(here("Outputs", "Final_TrajShape_Models", "Annual", "finalshape_annuals11_central.rds"))
am_trend_central<- am_trend_central %>% mutate(climate_zone = "central")


Sys.time(); step1(am_trend_central, "anisoEVI", "central"); Sys.time()
Sys.time(); step1(am_trend_southern, "anisoEVI", "southern"); Sys.time()
Sys.time(); step1(am_trend_eastern, "anisoEVI", "eastern"); Sys.time()

remove(am_trend_central, am_trend_eastern, am_trend_southern)

#EVI
am_trend_southern<- read_rds(here("Outputs", "Final_TrajShape_Models", "Annual_EVI", "finalshape_annuals11_southern.rds"))
am_trend_southern<- am_trend_southern %>% mutate(climate_zone = "southern")
am_trend_eastern<- read_rds(here("Outputs", "Final_TrajShape_Models", "Annual_EVI", "finalshape_annuals11_eastern.rds"))
am_trend_eastern<- am_trend_eastern %>% mutate(climate_zone = "eastern")
am_trend_central<- read_rds(here("Outputs", "Final_TrajShape_Models", "Annual_EVI", "finalshape_annuals11_central.rds"))
am_trend_central<- am_trend_central %>% mutate(climate_zone = "central")

Sys.time(); step1(am_trend_central, "EVI", "central"); Sys.time()
Sys.time(); step1(am_trend_southern, "EVI", "southern"); Sys.time()
Sys.time(); step1(am_trend_eastern, "EVI", "eastern"); Sys.time()

remove(am_trend_central, am_trend_eastern, am_trend_southern)

#NDVI
am_trend_southern<- read_rds(here("Outputs", "Final_TrajShape_Models", "Annual_NDVI", "finalshape_annuals11_southern.rds"))
am_trend_southern<- am_trend_southern %>% mutate(climate_zone = "southern")
am_trend_eastern<- read_rds(here("Outputs", "Final_TrajShape_Models", "Annual_NDVI", "finalshape_annuals11_eastern.rds"))
am_trend_eastern<- am_trend_eastern %>% mutate(climate_zone = "eastern")
am_trend_central<- read_rds(here("Outputs", "Final_TrajShape_Models", "Annual_NDVI", "finalshape_annuals11_central.rds"))
am_trend_central<- am_trend_central %>% mutate(climate_zone = "central")

Sys.time(); step1(am_trend_central, "NDVI", "central"); Sys.time()
Sys.time(); step1(am_trend_southern, "NDVI", "southern"); Sys.time()
Sys.time(); step1(am_trend_eastern, "NDVI", "eastern"); Sys.time()

remove(am_trend_central, am_trend_eastern, am_trend_southern)

##Note above resulting csvs were moved to server
##and the csv name was modifed manually to include _climatezone

#Match above dfs with stl dfs on server

step2<- function (sample_df, stl_df){
  cell_sample<- unique(sample_df$cell)
  x<- stl_df %>% filter(cell %in% cell_sample)
}

#Plotting time series for sampled pixels for each trajectory and each index from above
# csv files with stl time series were transferred back from server to appropriate folders manually
# and renamed files to sampled_stl_df so that I can run the below function

step3<- function (index_name, climate_zone){
  sampled_stl_df <- read.csv(here("Outputs", "TrajectoryPlotting", index_name, climate_zone, "sampled_stl_df.csv"))
  sampled_stl_df <- sampled_stl_df %>% rename("timestamp"  = "name")
  sampled_stl_df <- sampled_stl_df %>%  separate_wider_delim(cols = timestamp, delim = "_", names = c("Year","Month"))
  sampled_stl_df<-  sampled_stl_df %>% mutate(Time = paste0(Year, "/", Month, "/01")) %>% dplyr::select(-c(Year,Month))
  sampled_stl_df<- sampled_stl_df %>% mutate(Time = as.Date(Time))
  
  sampled_df<- read.csv(here("Outputs", "TrajectoryPlotting", index_name, climate_zone, "sampled_df.csv" ))
  sampled_df<- sampled_df %>%  mutate(value=case_when(pixelvalue==1~"Linear-decrease",
                                                      pixelvalue==2~"Linear-increase",
                                                      pixelvalue==4~"Step-decrease",
                                                      pixelvalue==5~"Step-increase",
                                                      pixelvalue==6~"Quad-decrease accelerated",
                                                      pixelvalue==7~"Quad-decrease deccelerated",
                                                      TRUE~"No trend"))
  
  unique_cells<- unique(sampled_df$cell)
  
  #plot
  plot_list<- list()
  for (i in 1:length(unique_cells)){
    unique_cell_stl<- sampled_stl_df %>% filter(cell == unique_cells[[i]])
    unique_cell_traj<- sampled_df %>% filter(cell == unique_cells[[i]])
    x_join <- left_join(unique_cell_stl, unique_cell_traj, by= c("cell" = "cell"))
    lp_fig<- x_join %>% 
      dplyr::select(c(value_int, trend)) %>%
      ggplot(aes(x = x_join$Time))+
      geom_line(aes(y = value_int), color = "darkred") + 
      geom_line(aes(y = trend), color="black", linetype="longdash")+
      theme_classic()+
      ylab("index") +
      xlab("Time")+
      ggtitle(unique(x_join$value.y))
    plot_list[[i]]<- lp_fig
  }
  
  x<- ggarrange(plotlist = plot_list)
  ggsave(here("Outputs", "TrajectoryPlotting", index_name, climate_zone, "plottedtimeseries.png"),x,
         dpi= 700, height = 70, width = 70, units= "cm")
}

Sys.time(); step3("NDVI", "central");Sys.time()
Sys.time(); step3("NDVI", "southern");Sys.time()
Sys.time(); step3("NDVI", "eastern");Sys.time()

Sys.time(); step3("EVI", "central");Sys.time()
Sys.time(); step3("EVI", "southern");Sys.time()
Sys.time(); step3("EVI", "eastern");Sys.time()

Sys.time(); step3("anisoEVI", "central");Sys.time()
Sys.time(); step3("anisoEVI", "southern");Sys.time()
Sys.time(); step3("anisoEVI", "eastern");Sys.time()

#Making map of above sampled points
step4<- function (index_name, climate_zone, d_trans_region_caponfirstletter, mappalette){
  sampled_df<- read.csv(here("Outputs", "TrajectoryPlotting", index_name, climate_zone, "sampled_df.csv" ))
  
  point_vector<- st_as_sf(terra::vect(sampled_df, geom=c("x", "y"), crs="epsg:4326"))
  point_vector<- point_vector %>%  mutate(value=case_when(pixelvalue==1~"Linear-decrease",
                                                      pixelvalue==2~"Linear-increase",
                                                      pixelvalue==4~"Step-decrease",
                                                      pixelvalue==5~"Step-increase",
                                                      pixelvalue==6~"Quad-decrease accelerated",
                                                      pixelvalue==7~"Quad-decrease deccelerated",
                                                      TRUE~"No trend"))
  
  background_shp<- d_trans %>% filter(region ==d_trans_region_caponfirstletter)
  
  map_points<-   
    tm_shape(background_shp)+ tm_fill(fill= "#CCCCCC")+
    tm_shape (point_vector)+
    tm_dots(size =0.7, fill = "value", fill.scale = tm_scale_categorical(values = mappalette)) 
  tmap_save(map_points, here("Outputs", "TrajectoryPlotting", index_name, climate_zone, "map_sampled_points.png"),
            height= 30, width= 30, units= "cm", dpi= 800)
}


mappalette1<- brewer.pal(n=7, "Set1")
mappalette1[[3]]<- "#333333"
mappalette1[[6]]<- "#984EA3"
mappalette1[[7]]<- "#FF7F00"
mappalette1[[4]]<- "#FFFF33"
mappalette1[[5]]<- "#A65628"
Sys.time(); step4("NDVI", "central", "Central", mappalette1);Sys.time() #because Quad trajectories exist

mappalette2<- brewer.pal(n=5, "Set1")
mappalette2[[3]]<- "#333333"
mappalette2[[4]]<- "#984EA3"
mappalette2[[5]]<- "#FF7F00"

Sys.time(); step4("NDVI", "southern", "Southern", mappalette2);Sys.time()
Sys.time(); step4("NDVI", "eastern", "Eastern", mappalette2);Sys.time()

Sys.time(); step4("EVI", "central", "Central", mappalette2);Sys.time()
Sys.time(); step4("EVI", "southern", "Southern", mappalette2);Sys.time()
Sys.time(); step4("EVI", "eastern", "Eastern", mappalette2);Sys.time()

Sys.time(); step4("anisoEVI", "central", "Central", mappalette2);Sys.time()
Sys.time(); step4("anisoEVI", "southern", "Southern", mappalette2);Sys.time()
Sys.time(); step4("anisoEVI", "eastern", "Eastern", mappalette2);Sys.time()



#ToDO 
#1) Complete stl central, southern and eastern bind_rows for NDVI- DONE
#2) Inner join by cell stl bind rows file to final results and check if works
#3) Take each raster produced at end of the for loop and exclude from trajmap df all pixels that have been masked out
#4) df from step 3, slice to extract 10 pixels per trajectory type
#5) Match to stl rows
#6) plot time series
#7) Make map of points selected 


##Extraction of time series for each index for different trajectories for same set of points
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


#Sampled points from NDVI index, per trajectory and in each climate zone
ndvi_sample_central<- read.csv(here("Outputs", "TrajectoryPlotting", "NDVI", "central", "sampled_df.csv"))
ndvi_sample_southern<- read.csv(here("Outputs", "TrajectoryPlotting", "NDVI", "southern", "sampled_df.csv"))
ndvi_sample_eastern<- read.csv(here("Outputs", "TrajectoryPlotting", "NDVI", "eastern", "sampled_df.csv"))

rename_ndvisample<- function (ndvi_sample_df){ #renaming to avoid confusion with other indices
  ndvi_sample_df<- ndvi_sample_df %>% dplyr::select(-c(X))
  names(ndvi_sample_df)[3]<- "NDVI_cellID"
  names(ndvi_sample_df)[5]<- "NDVI_traj"
  ndvi_sample_df
}

ndvi_sample_central<- rename_ndvisample(ndvi_sample_central)
ndvi_sample_southern<- rename_ndvisample(ndvi_sample_southern)
ndvi_sample_eastern<- rename_ndvisample(ndvi_sample_eastern)


NDVIsample_pointshp<- function(ndvi_sample_df, climate_zone){
  point_vector<- st_as_sf(terra::vect(ndvi_sample_df, geom=c("x", "y"), crs="epsg:4326"))
  st_write(point_vector, here("Outputs", "TrajectoryPlotting", "NDVI", paste0(climate_zone, "_NDVIsample.shp")))
} 

NDVIsample_pointshp(ndvi_sample_central, "central")
NDVIsample_pointshp(ndvi_sample_southern, "southern")
NDVIsample_pointshp(ndvi_sample_eastern, "eastern")
remove(ndvi_sample_central, ndvi_sample_eastern, ndvi_sample_southern, rename_ndvisample, NDVIsample_pointshp)

#Extracting cell index and traj from EVI and anisoEVI for above NDVI sampled points
step1<- function (final_results_df, index_name, climate_zone, ndvisample_pointshp){
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
  
  #raster of cell index
  cellindex_vector<- terra::vect(final_results_df, geom=c("x", "y"), crs="epsg:4326")
  Sys.time(); cellindex_raster<- terra::rasterize(cellindex_vector, evi, "cell", fun="max"); Sys.time()
  masked_cellindex_raster<- terra::mask(cellindex_raster, threshold_mapbiomas_mask_20)
  
  #sample NDVI points in above raster toi get cellID
  newindex_sample_cellID<- terra::extract(masked_cellindex_raster, ndvisample_pointshp, method ="simple", bind = TRUE)
  
  #for above cellID get trajectory from new index
  newindex_sample_cellID<- as.data.frame(newindex_sample_cellID)
  
    
  newindex_sample_celltraj<- left_join(newindex_sample_cellID,final_results_df, by= c("max"="cell"))
 
  
  names(newindex_sample_celltraj)[4]<-  paste0(index_name, "_ID")
  names(newindex_sample_celltraj)[8]<- paste0(index_name, "_tr")
  
  write.csv(newindex_sample_celltraj, here("Outputs", "TrajectoryPlotting", index_name, paste0("NDVI_", index_name, "_", climate_zone,  "_sampledpoints2.csv")))
  
}  

#Application of above function  
ndvi_sample_centralshp<- st_read(here("Outputs", "TrajectoryPlotting", "NDVI", "central_NDVIsample.shp"))
ndvi_sample_southernshp<- st_read(here("Outputs", "TrajectoryPlotting", "NDVI", "southern_NDVIsample.shp"))
ndvi_sample_easternshp<- st_read(here("Outputs", "TrajectoryPlotting", "NDVI", "eastern_NDVIsample.shp"))


#EVI
am_trend_southern<- read_rds(here("Outputs", "Final_TrajShape_Models", "Annual_EVI", "finalshape_annuals11_southern.rds"))
am_trend_southern<- am_trend_southern %>% mutate(climate_zone = "southern")
am_trend_eastern<- read_rds(here("Outputs", "Final_TrajShape_Models", "Annual_EVI", "finalshape_annuals11_eastern.rds"))
am_trend_eastern<- am_trend_eastern %>% mutate(climate_zone = "eastern")
am_trend_central<- read_rds(here("Outputs", "Final_TrajShape_Models", "Annual_EVI", "finalshape_annuals11_central.rds"))
am_trend_central<- am_trend_central %>% mutate(climate_zone = "central")

Sys.time(); step1(am_trend_central, "EVI", "central", ndvi_sample_centralshp); Sys.time()
Sys.time(); step1(am_trend_southern, "EVI", "southern", ndvi_sample_southernshp); Sys.time()
Sys.time(); step1(am_trend_eastern, "EVI", "eastern", ndvi_sample_easternshp); Sys.time()

remove(am_trend_central, am_trend_eastern, am_trend_southern)

#anisoEVI
am_trend_southern<- read_rds(here("Outputs", "Final_TrajShape_Models", "Annual", "finalshape_annuals11_southern.rds"))
am_trend_southern<- am_trend_southern %>% mutate(climate_zone = "southern")
am_trend_eastern<- read_rds(here("Outputs", "Final_TrajShape_Models", "Annual", "finalshape_annuals11_eastern.rds"))
am_trend_eastern<- am_trend_eastern %>% mutate(climate_zone = "eastern")
am_trend_central<- read_rds(here("Outputs", "Final_TrajShape_Models", "Annual", "finalshape_annuals11_central.rds"))
am_trend_central<- am_trend_central %>% mutate(climate_zone = "central")


Sys.time(); step1(am_trend_central, "anisoEVI", "central", ndvi_sample_centralshp); Sys.time()
Sys.time(); step1(am_trend_southern, "anisoEVI", "southern",  ndvi_sample_southernshp); Sys.time()
Sys.time(); step1(am_trend_eastern, "anisoEVI", "eastern", ndvi_sample_easternshp); Sys.time()

#Extracting time series from anisoEVI and EVI for above cellIDs
step2<- function (sample_df, stl_df){
  cell_sample<- unique(sample_df$cell)
  x<- stl_df %>% filter(cell %in% cell_sample)
}
 
#Extracted time series from stl files on server. Transferred back to this laptop the
#csv of the time series

#Preparing df for plotting
ndvi_timeseries_central <- read.csv(here("Outputs", "TrajectoryPlotting", "NDVI", "central", "sampled_stl_df.csv"))
ndvi_timeseries_southern <- read.csv(here("Outputs", "TrajectoryPlotting", "NDVI", "southern", "sampled_stl_df.csv"))
ndvi_timeseries_eastern <- read.csv(here("Outputs", "TrajectoryPlotting", "NDVI", "eastern", "sampled_stl_df.csv"))

step3a<- function (ndvi_timeseries, sampled_ndvi_evi_df){#connect ndvi time series to ndvi sampled points with EVI
  ndvi_timeseries <- ndvi_timeseries  %>% dplyr::select(c(cell, name, value_int, trend))
  ndvi_timeseries_join<- left_join(sampled_ndvi_evi_df, ndvi_timeseries, by= c("NDVI_ID" = "cell"))
  
  ndvi_timeseries_join<- ndvi_timeseries_join %>% rename("NDVI_value_int" = "value_int")
  ndvi_timeseries_join<- ndvi_timeseries_join %>% rename("NDVI_trend" = "trend")
  ndvi_timeseries_join
}

#EVI
sampled_points_ndvi_evi_central<- read.csv(here("Outputs", "TrajectoryPlotting", "EVI", "NDVI_EVI_central_sampledpoints2.csv"))
sampled_points_ndvi_evi_southern<- read.csv(here("Outputs", "TrajectoryPlotting", "EVI", "NDVI_EVI_southern_sampledpoints2.csv"))
sampled_points_ndvi_evi_eastern<- read.csv(here("Outputs", "TrajectoryPlotting", "EVI", "NDVI_EVI_eastern_sampledpoints2.csv"))

ndvitimeseries_evi_central<- step3a(ndvi_timeseries_central, sampled_points_ndvi_evi_central)
ndvitimeseries_evi_southern<- step3a(ndvi_timeseries_southern, sampled_points_ndvi_evi_southern)
ndvitimeseries_evi_eastern<- step3a(ndvi_timeseries_eastern, sampled_points_ndvi_evi_eastern)

remove(ndvi_timeseries_central, ndvi_timeseries_eastern, ndvi_timeseries_southern)

step3b<- function(evi_timeseries, ndvitimeseries_index, ndvitimseries_join_col ){#connect evi time series to previous dfs
  evi_cell_id<- unique(ndvitimeseries_index$ndvitimseries_join_col)
  
  joinedlist<- list()
  for (i in 1:length(evi_cell_id)){
    ndvitimeseries_subset<- ndvitimeseries_index %>% filter(ndvitimseries_join_col == evi_cell_id[[i]])
    evitimeseries_subset<- evi_timeseries %>% filter(cell == evi_cell_id[[i]])
    evitimeseries_subset<- evitimeseries_subset %>% dplyr::select(c(cell, value_int, trend))
    ndvits_evits<- bind_cols (ndvitimeseries_subset, evitimeseries_subset)
    ndvits_evits<- ndvits_evits %>% rename("EVI_value_int" = "value_int")
    ndvits_evits<- ndvits_evits %>% rename("EVI_trend" = "trend")
    joinedlist[[i]]<- ndvits_evits
    }
  ndvitimeseries_evitimeseries<- bind_rows(joinedlist)
  ndvitimeseries_evitimeseries
}

# EVI timeseries  
evits_central<- read.csv(here("Outputs", "TrajectoryPlotting", "EVI", "NDVI_EVI_central_sample_stl.csv"))
evits_southern<- read.csv(here("Outputs", "TrajectoryPlotting", "EVI", "NDVI_EVI_southern_sample_stl.csv"))
evits_eastern<- read.csv(here("Outputs", "TrajectoryPlotting", "EVI", "NDVI_EVI_eastern_sample_stl.csv"))

ndvi_evi_ts_central<- step3b(evits_central, ndvitimeseries_evi_central, EVI_ID)
ndvi_evi_ts_southern<- step3b(evits_southern, ndvitimeseries_evi_southern, EVI_ID)
ndvi_evi_ts_eastern<- step3b(evits_eastern, ndvitimeseries_evi_eastern, EVI_ID)

remove(evits_central, evits_eastern, evits_southern)
remove(sampled_points_ndvi_evi_central, sampled_points_ndvi_evi_eastern, sampled_points_ndvi_evi_southern)
remove(ndvitimeseries_evi_central, ndvitimeseries_evi_eastern, ndvitimeseries_evi_southern)

step3c<- function (anisoEVI_timeseries, ndvi_anisoEVI_sample){
  ndvi_anisoEVI_join<- left_join(anisoEVI_timeseries, ndvi_anisoEVI_sample, by= c("cell" = "anisoEVI_ID"))
  ndvi_anisoEVI_join<- ndvi_anisoEVI_join %>% dplyr::select(c(NDVI_ID, cell, x.x, y.x, name, value_int, trend, anisoEVI_tr))
  ndvi_anisoEVI_join<- ndvi_anisoEVI_join %>% rename("anisoEVI_ID" = "cell")
}

#anisoEVI
sampled_points_ndvi_anisoevi_central<- read.csv(here("Outputs", "TrajectoryPlotting", "anisoEVI", "NDVI_anisoEVI_central_sampledpoints2.csv"))
sampled_points_ndvi_anisoevi_southern<- read.csv(here("Outputs", "TrajectoryPlotting", "anisoEVI", "NDVI_anisoEVI_southern_sampledpoints2.csv"))
sampled_points_ndvi_anisoevi_eastern<- read.csv(here("Outputs", "TrajectoryPlotting", "anisoEVI", "NDVI_anisoEVI_eastern_sampledpoints2.csv"))

# anisoEVI timeseries  
anisoevits_central<- read.csv(here("Outputs", "TrajectoryPlotting", "anisoEVI", "NDVI_anisoEVI_central_sample_stl.csv"))
anisoevits_southern<- read.csv(here("Outputs", "TrajectoryPlotting", "anisoEVI", "NDVI_anisoEVI_southern_sample_stl.csv"))
anisoevits_eastern<- read.csv(here("Outputs", "TrajectoryPlotting", "anisoEVI", "NDVI_anisoEVI_eastern_sample_stl.csv"))

ndvi_anisoEVItimeseries_central <- step3c (anisoevits_central,sampled_points_ndvi_anisoevi_central)
ndvi_anisoEVItimeseries_southern <- step3c (anisoevits_southern,sampled_points_ndvi_anisoevi_southern)
ndvi_anisoEVItimeseries_eastern <- step3c (anisoevits_eastern,sampled_points_ndvi_anisoevi_eastern)

remove(sampled_points_ndvi_anisoevi_central, sampled_points_ndvi_anisoevi_eastern, sampled_points_ndvi_anisoevi_southern)
remove(anisoevits_central, anisoevits_eastern, anisoevits_southern)


step3d<- function (joined_evits_df, joined_anisoevits_df) { #final join
  df_list<- list()
  unique_NDVIid<- unique(joined_anisoevits_df$NDVI_ID)
  for (i in 1:length(unique_NDVIid)){
    evits<- joined_evits_df %>% filter(NDVI_ID== unique_NDVIid[[i]])
    anisoevits<- joined_anisoevits_df %>% filter(NDVI_ID == unique_NDVIid[[i]])
    full_join<- bind_cols(evits, anisoevits)
    full_join<- full_join %>% dplyr::select(NDVI_ID...2, clmt_zn, NDVI_tr, EVI_ID, x, y, EVI_tr, NDVI_value_int, NDVI_trend, EVI_value_int, EVI_trend, anisoEVI_ID, value_int, trend, anisoEVI_tr)
    full_join<- full_join %>% rename("NDVI_ID" = "NDVI_ID...2")
    full_join<- full_join %>% rename("anisoEVI_value_int" = "value_int")
    full_join<- full_join %>% rename("anisoEVI_trend" = "trend")
    df_list[[i]]<- full_join
  }
  final_df<- bind_rows(df_list)
  final_df
}

central_df<- step3d(ndvi_evi_ts_central, ndvi_anisoEVItimeseries_central)
eastern_df<- step3d(ndvi_evi_ts_eastern, ndvi_anisoEVItimeseries_eastern)
southern_df<- step3d(ndvi_evi_ts_southern, ndvi_anisoEVItimeseries_southern)

remove(ndvi_anisoEVItimeseries_central, ndvi_anisoEVItimeseries_eastern, ndvi_anisoEVItimeseries_southern)
remove(ndvi_evi_ts_central, ndvi_evi_ts_eastern, ndvi_evi_ts_southern)

#Plotting time series
final_df<- bind_rows(central_df, eastern_df, southern_df)
#write.csv(final_df, here("Outputs", "TrajectoryPlotting", "sampled_points_timeseries_df.csv"))

final_df2<- final_df %>% dplyr::select(c(NDVI_ID, NDVI_tr, EVI_tr, anisoEVI_tr, 
                                         NDVI_value_int, NDVI_trend, EVI_value_int, 
                                         EVI_trend, anisoEVI_value_int, anisoEVI_trend))
pivot_df<- pivot_longer(final_df2, 5:10)
pivot_df<- pivot_df %>% mutate(name = str_replace(name, "_value_int", "_value"))
pivot_df<- pivot_df %>% separate_wider_delim(cols = name, delim = "_", names = c("IndexName","value1"))

pivot_df2<- pivot_df %>% pivot_wider(names_from = value1, values_from = value)
pivot_df2<- pivot_df2 %>%  mutate(NDVI_tr = case_when(NDVI_tr ==1~"Linear-decrease",
                                                      NDVI_tr==2~"Linear-increase",
                                                      NDVI_tr==4~"Step-decrease",
                                                      NDVI_tr==5~"Step-increase",
                                                      NDVI_tr==6~"Quad-decrease accelerated",
                                                      NDVI_tr==7~"Quad-decrease deccelerated",
                                                                          TRUE~"No trend"))
pivot_df2<- pivot_df2 %>%  mutate(EVI_tr = case_when(EVI_tr ==1~"Linear-decrease",
                                                      EVI_tr==2~"Linear-increase",
                                                      EVI_tr==4~"Step-decrease",
                                                      EVI_tr==5~"Step-increase",
                                                      EVI_tr==6~"Quad-decrease accelerated",
                                                      EVI_tr==7~"Quad-decrease deccelerated",
                                                      TRUE~"No trend"))
pivot_df2<- pivot_df2 %>%  mutate(anisoEVI_tr = case_when(anisoEVI_tr ==1~"Linear-decrease",
                                                          anisoEVI_tr==2~"Linear-increase",
                                                          anisoEVI_tr==4~"Step-decrease",
                                                          anisoEVI_tr==5~"Step-increase",
                                                          anisoEVI_tr==6~"Quad-decrease accelerated",
                                                          anisoEVI_tr==7~"Quad-decrease deccelerated",
                                                      TRUE~"No trend"))



unique_ID<- unique(pivot_df2$NDVI_ID)
list_fig<- list()
for (i in 1:length(unique_ID)){
  
  df<- pivot_df2 %>% filter(NDVI_ID == unique_ID[[i]])
  df<- df%>% unnest (value, trend)
  
  
  lp_fig<- df %>% mutate(time= rep(seq(1:240), times=3)) %>%  
    ggplot(aes(x = time, color = IndexName))+
    geom_line(aes(y = value, alpha = 0.4)) + 
    geom_line(aes(y = trend), lwd= 1.5, linetype="longdash")+
    theme_classic()+  scale_color_manual(values = c("#648FFF","#DC267F", "#FE6100"))+
    ylab("index") + 
    xlab("Time") + scale_alpha(guide = 'none')+
    ggtitle (paste0("NDVI = ", unique(df$NDVI_tr), ";",
                    " EVI = ", unique(df$EVI_tr), ";",
                    " anisoEVI = ", unique(df$anisoEVI_tr)))
  
    list_fig[[i]]<- lp_fig
}

#Map making
centralNDVIsample<- st_read(here("Outputs","TrajectoryPlotting", "NDVI", "central_NDVIsample.shp"))
southernNDVIsample<- st_read(here("Outputs","TrajectoryPlotting", "NDVI", "southern_NDVIsample.shp"))
easternNDVIsample<- st_read(here("Outputs","TrajectoryPlotting", "NDVI", "eastern_NDVIsample.shp"))

NDVIsample<- bind_rows(centralNDVIsample, southernNDVIsample, easternNDVIsample)

mappalette<- c("red", "black")
map_list<- list()
for (i in 1:length(unique_ID)){
  all_points<- NDVIsample
  all_points<- all_points %>% mutate(colorcode = case_when(NDVI_ID==unique_ID[[i]] ~"1",
                                                           TRUE ~ "2"))

  map_points<-   
    tm_shape(cerrado)+ tm_fill(fill= "#CCCCCC")+
    tm_shape (all_points)+
    tm_dots(size =0.5, fill = "colorcode", 
            fill.scale = tm_scale_categorical(values = mappalette)) +
    tm_layout(legend.show = FALSE) 
  map_list[[i]]<- map_points
  
}

for ( i in 1:length(unique_ID)){
  grobobject<- tmap_grob(map_list[[i]])
  fig<- ggarrange(grobobject, list_fig[[i]], widths= c(1,2))
  ggsave(here("Outputs", "TrajectoryPlotting", "allindices_timeseries_plots_sampled_points",
              paste0("plot_", unique_ID[[i]], ".png")), fig,
         dpi= 700, height = 10, width = 30, units= "cm")
  print (i)
}


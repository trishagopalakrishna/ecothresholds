library(terra)
library(sf)
library(tmap)
library(RColorBrewer)
library(here)
library(tidyverse)


#Main figure- NDVI monthly unstable
annual_ndvi_5km_unstablearea <- rast(here("Outputs", "TrendsResults", "results_rasters" , "annual", "unstablearea_5km", "annualswin11_ndvi_unstablearea.tif"))
monthly_ndvi_5km_unstablearea <- rast(here("Outputs", "TrendsResults", "results_rasters" , "monthly", "unstablearea_5km" ,"monthlyswin11_ndvi_unstablearea.tif"))

annual_evi_5km_unstablearea <- rast(here("Outputs", "TrendsResults", "results_rasters" , "annual" ,"unstablearea_5km", "annualswin11_evi_unstablearea.tif"))
monthly_evi_5km_unstablearea <- rast(here("Outputs", "TrendsResults", "results_rasters" , "monthly" ,"unstablearea_5km","monthlyswin11_evi_unstablearea.tif"))


annual_aniso_5km_unstablearea <- rast(here("Outputs", "TrendsResults", "results_rasters" , "annual", "unstablearea_5km","annualswin11_aniso_unstablearea.tif"))
monthly_aniso_5km_unstablearea <- rast(here("Outputs", "TrendsResults", "results_rasters" , "monthly" , "unstablearea_5km","monthlyswin11_aniso_unstablearea.tif"))


d_trans<- st_read(here ("Data", "Cattelanetal_Clustering", "cerrado_climate_zones.shp"))
d_trans<- st_transform(d_trans, crs = 4326)
d_trans<- d_trans %>% mutate(area_km= terra::expanse(vect(d_trans), unit="km"))
cerrado<- d_trans %>% st_union()

static_map_function<- function (raster,index_name){
  s_map<-
    tm_shape(cerrado)+ 
    tm_borders()+
    tm_shape (raster)+
    tm_layout(frame= F)+
    tm_raster(col.scale = tm_scale_continuous(values ="rd_pu"),
              col.legend = tm_legend(position = c("left", "top"), title ="Proportion of 5x5km pixel",
                                     orientation = "landscape", frame.lwd = 0.2, width = 12, text.size = 3, title.size =3)) +
    tm_title (index_name, size = 5, position = tm_pos_out("left", "top")) 
  s_map
}

ndvi_monthly<- static_map_function(monthly_ndvi_5km_unstablearea, "NDVI")
tmap_save(ndvi_monthly, filename = here("Outputs", "ATBC2025","NDVI_monthly_unstablearea2.jpg"), 
          height = 80, width = 80, units = "cm")



# Barplot _ NDVI monthly trajectories barplot
monthly_input_file_path <- here("Outputs", "TrendsResults", "results_rasters", "monthly")
monthly_file_list<- list.files(path = paste0(monthly_input_file_path, "/"), pattern = paste0("*","_monthly_swin11","*"), all.files = T, full.names = T)
monthly_file_list <- gtools::mixedsort(monthly_file_list)

monthly_rds_list<- lapply(monthly_file_list, rast)
names(monthly_rds_list[[2]])<- "evi_swin11_monthly"
names(monthly_rds_list[[3]])<- "ndvi_swin11_monthly"
names(monthly_rds_list[[1]])<- "aniso_swin11_monthly"

df_prep_function<- function (trendresults_raster, index_name){
  x_df<- terra::as.data.frame(trendresults_raster, cell = TRUE, xy= TRUE)
  x_df<- x_df %>% mutate(index = index_name)
  names(x_df)[4]<- "value"
  percentage_df <- x_df %>% group_by(value, index) %>% summarise(count= n()) %>% mutate(percentage= (count/nrow(x_df))*100) 
  percentage_df
}

Sys.time(); monthly_ndvi_trend_results <- df_prep_function(monthly_rds_list[[3]], "NDVI (monthly)"); Sys.time()

plot_df <- bind_rows(monthly_ndvi_trend_results)

plot_df <- plot_df %>% mutate(value = case_when(value== 1~ "Linear decrease",
                                                value==2~ "Linear increase",
                                                value==4~ "Step decrease",
                                                value==5~"Step increase",
                                                value==6~"Quadratic decrease (accelerated)",
                                                value==7~"Quadratic decrease (decelerated)",
                                                value==8~"Quadratic increase (accelerated)",
                                                value==9~"Quadratic increase (decelerated)",
                                                TRUE~"No trend"))
plot_df <- plot_df %>%  mutate(my_alpha = ifelse(value == "Step increase" | value == "No trend", 1, 0.5))

Set1palette<- c("#44AA99", "#332288", "slateblue1", "#CC6677", "#AA4499", "#88CCEE", "#DDCC77" , "#882255",  "sienna3")

trendresults_barplot<- plot_df %>% 
  ggplot(aes(value, percentage, fill = value, group=index, alpha=my_alpha)) +
  geom_col(position = position_dodge(), color = "black") +
  theme_classic() +
  ylab ("% of native savannas") +
  theme(strip.placement  = "outside",
        panel.spacing    = unit(0, "points"),
        strip.background = element_blank(),
        strip.text       = element_text(size = 8, angle = 90),
        axis.text.x = element_text(size = 8, angle = 90, hjust = .5, vjust = .5),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = Set1palette) + theme(legend.position="none") 
ggsave(here("Outputs", "ATBC2025","NDVI_monthly_trajectories_barplot.jpg"),
       trendresults_barplot,dpi = 700, height = 10, width=10, units = "cm")


#Animations of monthly unstable area for all indices
index_monthly_unstable<- c(monthly_ndvi_5km_unstablearea, monthly_evi_5km_unstablearea, monthly_aniso_5km_unstablearea)
names(index_monthly_unstable)<- c("NDVI", "EVI", "anisotropic EVI")
animation_map<-
    tm_shape(cerrado)+ tm_borders()+
  tm_layout(frame= F, title.position = tm_pos_in("left", "top"), panel.label.frame = F, panel.label.bg.color = "white") +
    tm_shape (index_monthly_unstable)+
    tm_raster(col.scale = tm_scale_continuous(values ="rd_pu"),
              col.legend = tm_legend(position = c("left", "top"), title ="Proportion of 5x5km pixel",
                                     orientation = "landscape", frame.lwd = 0.2, width = 12, text.size = 0.5, title.size =0.5)) +
    tm_facets(nrow=1, ncol=1, free.coords = TRUE) 
tmap_animation(animation_map, filename = here("Outputs", "ATBC2025","trial_animation.gif"),
               delay = 50, loop = TRUE, restart.delay = 500,
               width = 800, height = 800, dpi = 150)

#Example NDVI time series of points that are step increase in monthly, but no trend in annual
ndvi_evi_ts_complete_tr <- read.csv(here("Outputs", "TrajectoryPlotting", "ndvi_evi_ts_trajectories.csv"))

eg1 <- ndvi_evi_ts_complete_tr %>% filter(NDVI_ID == 1293359) %>%
  dplyr::select(c(NDVI_original, NDVI_trend, time))
pivot_eg1 <- eg1 %>% pivot_longer(1:2)
pivot_eg1 <- pivot_eg1 %>% mutate(name = case_when(name == "NDVI_original"~ "Original",
                                                   TRUE ~ "De-seasoned")) 

lp_fig1<- pivot_eg1  %>%  
  ggplot(aes(x = as.factor(time)))+
  geom_line(aes(y = value, group= name, linetype = name))+
  scale_linetype_manual(values=c("solid", "dotted")) +
  theme_classic() +
  theme (axis.text.x = element_blank(),
         axis.ticks = element_blank(), 
         axis.title.x = element_blank(),
         legend.title = element_blank(),
         axis.title.y = element_blank(), legend.position = "none") 
ggsave(here("Outputs", "ATBC2025","NDVI_timeseries_lineplot.jpg"),
       lp_fig1,dpi = 700, height = 8, width=8, units = "cm")

#plot studies lollipop plot for woody encroachment

stevens_si<- read_csv(here("Data","Stevensetal_woodyencroachment_metaanalyses_data", "si_plotinformation.csv"))
names(stevens_si)[14]<-"y"
names(stevens_si)[15]<-"x"
stevens_si_removeNA<- stevens_si %>% filter(!is.na(x))
stevens_si_removeNA<- stevens_si_removeNA %>% filter(Country!="Venezuela")

#Duplicate plots resulted in multiple papers eg Durigan et al 1987 and 2006 are from the same plots
#So I "merge" rows and take the mean of the "annual change" that was calculated in Stevens et al 
x<-stevens_si_removeNA %>% 
  group_by(x,y) %>%
  summarise(mean_annual_rate= mean(Annual_change)) %>%
  ungroup()

x_join<- stevens_si_removeNA %>% left_join(x)
x_join2<- x_join %>% filter(!duplicated(mean_annual_rate))

stevens_si_vector<- sf::st_as_sf(x_join2, coords = c("x", "y"), crs="epsg:4326")
stevens_si_vector<- stevens_si_vector %>% mutate(uniqueID= 1:nrow(stevens_si_vector))


input_file_path <- here("Outputs", "TrendsResults", "results_rasters", "monthly")
file_list<- list.files(path = paste0(input_file_path, "/"), pattern = paste0("*","monthly_swin11","*"), all.files = T, full.names = T)
file_list <- gtools::mixedsort(file_list)
monthly_rds_list<- lapply(file_list, rast)

trajectory_segregate <- function(x) {
  x_split<- terra::segregate(x, other=999)
  x_split_crop<- terra::crop(x_split, vect(d_trans))
  x_split_mask<- terra::mask(x_split_crop, vect(d_trans))
  x_split_mask
}
monthly_NDVItrajectories_split<- trajectory_segregate(monthly_rds_list[[1]])

distance_function<- function(split_trajectory_raster){
  x_dist<- terra::gridDist(split_trajectory_raster, target= 1)
}

calculate_distances<- function (split_raster){
  list_distance_raster<- list()
  for (i in 1:nlyr(split_raster)){
    x_distance<- distance_function(split_raster[[i]])
    list_distance_raster[[i]]<- x_distance
  }
  distance_raster<- rast(list_distance_raster)
  names(distance_raster)<- 1:nlyr(split_raster)
  distance_raster
}

Sys.time(); monthly_NDVI_distances<- calculate_distances(monthly_NDVItrajectories_split); Sys.time() 

shortestdistance_trajectory<- function (distance_raster){
  distance_df<- terra::extract(distance_raster, stevens_si_vector)
  pivot_dist_df<- distance_df %>% pivot_longer(-1)
  pivot_dist_df2<-pivot_dist_df  %>% 
    group_by(as.factor(ID)) %>% 
    summarise(min(value, na.rm= TRUE), name[which.min(value)])
  names(pivot_dist_df2)<- c("uniqueID", "min_dist_m", "respective_trajshape")
  pivot_dist_df2
}
Sys.time(); monthly_NDVI_closest_trajectory_df<- shortestdistance_trajectory(monthly_NDVI_distances); Sys.time()

convert_NA_to_999_function <- function (split_trajectory_raster){
  classified_raster <- terra::classify(split_trajectory_raster, cbind (NA, 999))
  crop_classified_raster <- terra::crop (classified_raster, vect(cerrado))
  mask_classified_raster <- terra::mask (crop_classified_raster, vect(cerrado))
  mask_classified_raster
}

converted_monthly_NDVItrajectories_split <- convert_NA_to_999_function(monthly_NDVItrajectories_split)
Sys.time(); converted_monthly_NDVI_distances<- calculate_distances(converted_monthly_NDVItrajectories_split); Sys.time() 
Sys.time(); converted_monthly_NDVI_closest_trajectory_df<- shortestdistance_trajectory(converted_monthly_NDVI_distances); Sys.time()

#Studies/field sites that are within my study region
df_prep<- function (closest_distance_df, index_name){
  x<- closest_distance_df %>% 
    dplyr::select(-c(min_dist_m))
  names(x)[2]<- index_name
  x
}

monthly_closest_ndvi<- df_prep(monthly_NDVI_closest_trajectory_df, "monthly_NDVI") 
pivot_closest_df <- pivot_longer (monthly_closest_ndvi, 2)
pivot_closest_df <- pivot_closest_df %>% separate_wider_delim(cols = name, delim = "_", names = c("monthly","IndexName"))
pivot_closest_df <- pivot_closest_df %>% mutate(value = case_when(value== 1~ "Linear decrease",
                                                                  value==2~ "Linear increase",
                                                                  value==4~ "Step decrease",
                                                                  value==5~"Step increase",
                                                                  value==6~"Quadratic decrease (accelerated)",
                                                                  value==7~"Quadratic decrease (decelerated)",
                                                                  value==8~"Quadratic increase (accelerated)",
                                                                  value==9~"Quadratic increase (decelerated)",
                                                                  TRUE~"No trend"))  
pivot_closest_df <- pivot_closest_df %>% dplyr::select(-c(monthly, IndexName))
#Studies/field sites that are outside the study area of this project
converted_monthly_NDVI_closest_trajectory_df
names(converted_monthly_NDVI_closest_trajectory_df)<- c("uniqueID", "monthly_NDVI_mindist", "monthly_NDVI")
converted_monthly_NDVI_closest_trajectory_df <- converted_monthly_NDVI_closest_trajectory_df %>% dplyr::select(-c(monthly_NDVI_mindist))
pivot_closest_df2 <- pivot_longer (converted_monthly_NDVI_closest_trajectory_df, 2)
pivot_closest_df2 <- pivot_closest_df2 %>% mutate(value = case_when(value== 1~ "Linear decrease",
                                                                  value==2~ "Linear increase",
                                                                  value==4~ "Step decrease",
                                                                  value==5~"Step increase",
                                                                  value==6~"Quadratic decrease (accelerated)",
                                                                  value==7~"Quadratic decrease (decelerated)",
                                                                  value==8~"Quadratic increase (accelerated)",
                                                                  value==9~"Quadratic increase (decelerated)",
                                                                  TRUE~"No trend")) 
pivot_closest_df2 <- pivot_closest_df2 %>% dplyr::select(-c(name))
plot_df <- rbind(pivot_closest_df, pivot_closest_df2)
plot_df2 <- plot_df %>% group_by(value) %>% summarise(n= n())


x <- ggplot(plot_df2, aes(x=value, y=n)) +
  geom_segment( aes(x=value, xend=value, y=0, yend=n), color="grey") +
  geom_point( color="orange", size=4) +
  theme_classic() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("") +
  ylab("# woody encroachment studies \n (Stevens et al., 2017)")
ggsave(here("Outputs", "ATBC2025","woodyencroachment_plot.jpg"),
       x,dpi = 700, height = 8, width=8, units = "cm")


# CVNP year of step increase- 2017
#ran all chunks in 1_stepchangepixel_analyses.Rmd in cvnp_analyses folder
#in last chunk I stopped just before the barplot lines ie line 210
barplot_df

ndvi_step <- barplot_df %>% 
  filter(StepType == "increase" & name == "NDVI(monthly)") %>% 
  dplyr::select(c(value, percentage_totalstep)) 

#ran lines 230 -266 from above script
barplot_df

ndvi_step2 <- barplot_df %>% 
  filter(StepType == "increase" & name == "NDVI(monthly)") %>% ungroup() %>%
  dplyr::select(c(month,name, percentage_totalstep))
names(ndvi_step2)[1]<- "value"
ndvi_step2<- ndvi_step2 %>% mutate(name ="ND")

ndvi_step$value<- as.factor(ndvi_step$value)
plot_df <- rbind(ndvi_step, ndvi_step2)

cols <- plot_df %>% 
  mutate(col = ifelse(value == 2018 | (value == "Apr" | value== "May" | value =="Jun" | value =="Jul"), "#CC6600", "lightgrey")) %>%
  ungroup () %>%
  dplyr::select(-c(name, value, percentage_totalstep)) 
cols <- as.vector(cols['col'])

y<-ggplot(data= plot_df, aes(x= name, y= percentage_totalstep, fill= as.factor(value))) +
  geom_bar(stat="identity", width= 0.2) + 
  scale_fill_manual(values = cols$col) +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none",
        axis.text.y= element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
ggsave(here("Outputs", "ATBC2025","cvnp_yearstep.jpg"),
       y,dpi = 700, height = 4, width=4, units = "cm")

# CVNP veg proportions
#In cvnp_analyses>3_trajectory_vegproportions, I ran the first and second chunks
data_only_ndvimonthly <- c(coarse_lev1_phy, cvnp_monthly_ndvi)
df_data <- terra::as.data.frame(data_only_ndvimonthly, cells = TRUE, xy = TRUE)
df_data <- df_data %>% rename("cell" ="cell", "x"="x", "y"="y",
                              "percentage_grass"="1", "percentage_savanna"="2",
                              "percentage_woodland"= "3", "percentage_anthropic"="4",
                              "monthly_NDVI"="ndvi_swin11_monthly")
summary(df_data)

df_data_noNA <- df_data %>% filter(!is.na(monthly_NDVI)) #NA values are around the boundary and a couple of water pixels
df_data_noNA  <- df_data_noNA  %>% mutate(TrajType = case_when(monthly_NDVI == 2~" Linear Increase",
                                                               monthly_NDVI == 3 ~"No trend",
                                                               monthly_NDVI == 4 ~"Step Decrease",
                                                               monthly_NDVI == 5 ~ "Step Increase",
                                                               monthly_NDVI == 8 ~ "Quadratic Increase (acc)",
                                                               monthly_NDVI == 9 ~ "Quadratic Increase (dec)"))
df_data_noNA <- df_data_noNA %>%  mutate(my_alpha = ifelse(TrajType == "Step Increase" | TrajType == "No trend", 1, 0.5))

plot<- ggplot(df_data_noNA, aes(x = reorder(TrajType, percentage_savanna), y = percentage_savanna, fill = TrajType, alpha = my_alpha)) + 
  geom_boxplot() + ylab("Proportion covered by savanna physionomies")+
  theme_classic(base_size = 14) +
  scale_fill_manual(values=c("#44AA99", "slateblue1", "#88CCEE", "#DDCC77" , "#882255",  "sienna3")) +
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
    axis.text.x =element_blank(),
    axis.ticks.x =element_blank())
ggsave(here("Outputs", "ATBC2025","cvnp_trajectory_savanna.jpg"),
       plot,dpi = 700, height = 15, width= 15, units = "cm")

# CVNP heterogeniety 
#In cvnp_analyses>4_trajectory_heterogeniety, I ran the first and second chunk (only till the map making lines ie line 81)
data_stack <- c(cvnp_monthly_ndvi, heterogen)
data_df <- terra::as.data.frame(data_stack)
data_df_noNA <- data_df %>% drop_na()
names(data_df_noNA) <- c("monthly_NDVI", "heterogeneity")

data_df_noNA <- data_df_noNA %>% mutate(TrajType = case_when(monthly_NDVI == 2~" Linear Increase",
                                                             monthly_NDVI == 3 ~"No trend",
                                                             monthly_NDVI == 4 ~"Step Decrease",
                                                             monthly_NDVI == 5 ~ "Step Increase",
                                                             monthly_NDVI == 8 ~ "Quadratic Increase (acc)",
                                                             monthly_NDVI == 9 ~ "Quadratic Increase (dec)"))


data_df_noNA  <- data_df_noNA %>%  mutate(my_alpha = ifelse(TrajType == "Step Increase" | TrajType == "No trend", 1, 0.5))

plot2<- ggplot(data_df_noNA, aes(x = reorder(TrajType, heterogeneity, FUN = median), y = heterogeneity, fill = TrajType, alpha = my_alpha)) + 
  geom_boxplot() + ylab("native vegetation heterogeneity") +
  theme_classic(base_size = 14) +
  scale_fill_manual(values=c("#44AA99", "slateblue1", "#88CCEE", "#DDCC77" , "#882255",  "sienna3")) +
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x =element_blank(),
        axis.ticks.x =element_blank())
ggsave(here("Outputs", "ATBC2025","cvnp_trajectories_heterogeneity.jpg"),
       plot2,dpi = 700, height = 15, width= 15, units = "cm")

#cvnp burned area
#In cvnp_analyses>5_trajectory_burnedarea, I ran the first and second chunks and in the third chunk line 55 only
data_stack <- c(cvnp_monthly_ndvi, mean_cvnp_burnedarea)
data_df <- terra::as.data.frame(data_stack)
data_df_noNA <- data_df %>% drop_na()
names(data_df_noNA) <- c("monthly_NDVI", "average_burnedarea")

data_df_noNA <- data_df_noNA%>% mutate(TrajType = case_when(monthly_NDVI == 2~" Linear Increase",
                                                            monthly_NDVI == 3 ~"No trend",
                                                            monthly_NDVI == 4 ~"Step Decrease",
                                                            monthly_NDVI == 5 ~ "Step Increase",
                                                            monthly_NDVI == 8 ~ "Quadratic Increase (acc)",
                                                            monthly_NDVI == 9 ~ "Quadratic Increase (dec)"))

data_df_noNA <- data_df_noNA %>%  mutate(my_alpha = ifelse(TrajType == "Step Increase" | TrajType == "No trend", 1, 0.5))

plot3<- ggplot(data_df_noNA, aes(x = reorder(TrajType, average_burnedarea, FUN = median), y = average_burnedarea, fill = TrajType, alpha = my_alpha)) + 
  geom_boxplot() + ylab("average %burnedarea") + 
  theme_classic(base_size = 14) +
  scale_fill_manual(values=c("#44AA99", "slateblue1", "#88CCEE", "#DDCC77" , "#882255",  "sienna3")) +
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x =element_blank(),
        axis.ticks.x =element_blank())
ggsave(here("Outputs", "ATBC2025","cvnp_trajectory_averageburnedarea.jpg"),
       plot3,dpi = 700, height = 15, width= 15, units = "cm")

library(tidyverse)
library(terra)
library(foreach)
library(doParallel)
library(mgcv)
library(gtools)

######################## Data input
index_file_path <- "/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/GAM_data/"

#response- binary trajectory (all 9) from NDVI swin 11 stl setting
monthly_input_file_path <- paste0(index_file_path, "TrendsResults", "results_rasters", "monthly", "binary_count_5km")
monthly_file_list<- list.files(path = monthly_input_file_path, pattern = paste0("*","ndvi_","*"), all.files = T, full.names = T)
monthly_file_list <- gtools::mixedsort(monthly_file_list)
monthly_ndvi_trajectories<- lapply(monthly_file_list, rast)
remove(monthly_file_list, monthly_input_file_path)

#static drivers
anthropic_dist<- rast((index_file_path, "OtherVariables", "AnthropicDist", "5km_meandist_20thresholdsanthropicpixel.tif"))
hand<- rast(paste0(index_file_path, "OtherVariables","HAND", "cerrado_hand_5km.tif"))

#trend variables
trend_burnedarea<- rast(paste0(index_file_path, "OtherVariables", "Fire", "theilsen_burnedarea_5km.tif"))
trend_re<- rast(paste0(index_file_path, "OtherVariables", "Climate", "theilsen_re.tif"))
trend_at<- rast(paste0(index_file_path, "OtherVariables", "Climate", "theilsen_at_5km.tif"))
trend_spei<-rast(paste0(index_file_path, "OtherVariables", "Climate", "theilsen_spei.tif"))
trend_heatwave<- rast(paste0(index_file_path, "OtherVariables", "Climate", "theilsen_heatwaves.tif"))

#mean variables
mean_burnedarea<- rast(paste0(index_file_path, "OtherVariables", "Fire", "mean_burnedarea_5km_2002_2021.tif"))
mean_re <- rast(paste0(index_file_path, "OtherVariables", "Climate","mean_annual_relative_entropy_1981_2021.tif") )
mean_at <- rast(paste0(index_file_path, "OtherVariables", "Climate", "mat_1981_2021_5km.tif"))
mean_spei <- rast( paste0(index_file_path, "OtherVariables", "Climate", "mean_annualspei.tif")) 
mean_heatwave <- rast(paste0(index_file_path, "OtherVariables", "Climate", "mean_annual_heatwaves_2002_2021.tif"))

#trend in % area of veg formations
trend_forestpercentage <- rast(paste0(index_file_path, "OtherVariables", "Formations_Heterogeniety", "theilsen_forestpercentage_5km.tif"))
trend_savannapercentage <- rast(paste0(index_file_path, "OtherVariables", "Formations_Heterogeniety", "theilsen_savannapercentage_5km.tif"))
trend_grasspercentage <- rast(paste0(index_file_path, "OtherVariables", "Formations_Heterogeniety", "theilsen_grasspercentage_5km.tif"))

#mean %area of veg formation
mean_savannapercentage <- rast(paste0(index_file_path, "OtherVariables", "Formations_Heterogeniety", "mean_savanna_percentage_5km_2002_2021.tif"))
mean_grasspercentage <- rast(paste0(index_file_path, "OtherVariables", "Formations_Heterogeniety", "mean_grass_percentage_5km_2002_2021.tif"))
mean_forestpercentage <- rast(paste0(index_file_path, "OtherVariables", "Formations_Heterogeniety", "mean_forest_percentage_5km_2002_2021.tif"))


######################## Data structuring
crop_mask_function <- function(raster_to_be_cropmask){
  x_mask <- terra::mask(raster_to_be_cropmask, monthly_ndvi_trajectories[[1]][[1]])
  x_mask
}

cm_anthropic_dist <- crop_mask_function(anthropic_dist)
cm_hand <- crop_mask_function(hand)
cm_trend_burnedarea <- crop_mask_function(trend_burnedarea)
cm_trend_re <- crop_mask_function(trend_at)
cm_trend_at <- crop_mask_function(trend_at)
cm_trend_spei <- crop_mask_function(trend_spei)
cm_trend_heatwave <- crop_mask_function(trend_heatwave)
cm_mean_burnedarea <- crop_mask_function(mean_burnedarea)
cm_mean_re <- crop_mask_function(mean_re)
cm_mean_at <- crop_mask_function(mean_at)
cm_mean_spei <- crop_mask_function(mean_spei)
cm_mean_heatwave <- crop_mask_function(mean_heatwave)
cm_trend_forest <- crop_mask_function(trend_forestpercentage)
cm_trend_savanna <- crop_mask_function(trend_savannapercentage)
cm_trend_grass <- crop_mask_function(trend_grasspercentage)
cm_mean_forest <- crop_mask_function(mean_forestpercentage)
cm_mean_grass <- crop_mask_function(mean_grasspercentage)
cm_mean_savanna <- crop_mask_function(mean_savannapercentage)

constantdrivers <- c(cm_anthropic_dist, cm_hand, 
                     cm_trend_burnedarea, cm_trend_re, cm_trend_at, cm_trend_spei, cm_trend_heatwave, 
                     cm_mean_burnedarea, cm_mean_re, cm_mean_at, cm_mean_spei, cm_mean_heatwave)

df_prep <- function (trend_veg_formation, mean_veg_formation, colname_trend, colname_mean){
  x_stack <- c(monthly_ndvi_trajectories[[1]], constantdrivers, trend_veg_formation, mean_veg_formation)
  names(x_stack) <- c("ndvi_monthly_lin_dec", 
                      "ndvi_monthly_lin_inc",
                      "ndvi_monthly_stable",
                      "ndvi_monthly_step_dec",
                      "ndvi_monthly_step_inc",
                      "ndvi_monthly_quad_dec_acc",
                      "ndvi_monthly_quad_dec_decc",
                      "ndvi_monthly_quad_inc_acc",
                      "ndvi_monthly_quad_inc_dec",
                      "anthropic_dist", 
                      "hand", 
                      "trend_burnedarea",
                      "trend_re", 
                      "trend_annualtemp", 
                      "trend_spei",
                      "trend_heatwave", 
                      "mean_burnedarea",
                      "mean_re",
                      "mean_at",
                      "mean_spei",
                      "mean_heatwave",
                      colname_trend,
                      colname_mean
  )
  analyses_df<- terra::as.data.frame(x_stack, xy = TRUE)
  analyses_df
}

grass_df <- df_prep(trend_grasspercentage, mean_grasspercentage, "trend_grasspercentage", "mean_grasspercentage")
forest_df <- df_prep(trend_grasspercentage, mean_grasspercentage, "trend_forestpercentage", "mean_forestpercentage")
savanna_df <- df_prep(trend_grasspercentage, mean_grasspercentage, "trend_savannapercentage", "mean_savannapercentage")

remove(anthropic_dist, hand, trend_burnedarea, trend_re, trend_at, trend_spei, trend_heatwave,
       mean_burnedarea, mean_re, mean_at, mean_spei, mean_heatwave,
       trend_forestpercentage, trend_savannapercentage, trend_grasspercentage, 
       mean_savannapercentage, mean_grasspercentage, mean_forestpercentage)
remove(cm_anthropic_dist, cm_hand, cm_trend_burnedarea, cm_trend_re, cm_trend_at, cm_trend_spei, cm_trend_heatwave,
       cm_mean_burnedarea, cm_mean_re, cm_mean_at, cm_mean_spei, cm_mean_heatwave,
       cm_trend_forest, cm_trend_savanna, cm_trend_grass, 
       cm_mean_savanna, cm_mean_grass, cm_mean_forest)


remove_NA <- function (df_analyses){
  analyses_df_noNA <- df_analyses %>% drop_na() 
  analyses_df_noNA
}
grass_df_NA <- remove_NA(grass_df)
forest_df_NA <-remove_NA(forest_df)
savanna_df_NA <- remove_NA(savanna_df)

remove_HAND_outliers <- function (df_analyses_NA){
  no_outliers_sample_df <- df_analyses_NA %>% 
    mutate (IQR = IQR(hand),
            O_upper = quantile (hand, probs= c(0.75), na.rm= FALSE)+ 1.5*IQR,
            O_lower = quantile(hand, probs= c(0.25), na.rm= FALSE) - 1.5*IQR
    ) %>% 
    filter(O_lower <= hand & hand <= O_upper) %>%
    dplyr::select(-c(IQR, O_upper, O_lower))
  no_outliers_sample_df
}
grass_df_noNA_nooutliers <- remove_HAND_outliers(grass_df_NA)
forest_df_noNA_nooutliers <- remove_HAND_outliers(forest_df_NA)
savanna_df_noNA_nooutliers <- remove_HAND_outliers(savanna_df_NA)

remove(grass_df_NA, savanna_df_NA, forest_df_NA)
remove(crop_mask_function, df_prep, remove_HAND_outliers, remove_NA)

######################## GAM
summary(grass_df_noNA_nooutliers)
summary(forest_df_noNA_nooutliers)
summary(savanna_df_noNA_nooutliers)

Sys.time(); x_model <- mgcv::bam( cbind(ndvi_monthly_step_inc, 25- ndvi_monthly_step_inc) ~ s(anthropic_dist, bs= "ts") +
                                          s(hand, bs= "ts") +
                                          s(trend_re, bs= "ts") +
                                          s(trend_annualtemp, bs= "ts") +
                                          s(trend_heatwave, bs ="ts") +
                                          s(trend_spei, bs= "ts") +
                                          s(trend_burnedarea, bs = "ts") +
                                          s(trend_forestpercentage, bs= "ts") +
                                          
                                          s(mean_re, bs= "ts") +
                                          s(mean_at, bs= "ts") +
                                          s(mean_heatwave,  bs= "ts") +
                                          s(mean_spei, bs= "ts") +
                                          s(mean_burnedarea, bs= "ts") +
                                          s(mean_forestpercentage, bs= "ts") +
                                          
                                          s(x,y, bs= "ts") +
                                          
                                          ti(mean_at, x,y, d=c(1,2), bs=c("ts", "ts"), k=c(5,12)) +
                                          ti(trend_annualtemp, x,y, d=c(1,2), bs=c("ts", "ts"), k=c(5,12)) +
                                          ti(trend_annualtemp, trend_re, bs=c("ts", "ts"), k=c(5,5)) +
                                          ti(trend_heatwave, trend_spei, bs=c("ts", "ts"), k=c(5,5)) +
                                          ti(anthropic_dist, trend_burnedarea, bs=c("ts", "ts"), k=c(5,5)) +
                                          ti(anthropic_dist, trend_forestpercentage, bs=c("ts", "ts"), k=c(5,5)) + 
                                          ti(trend_burnedarea, trend_annualtemp, trend_re, bs=c("ts", "ts", "ts"), k=c(5,5,5)) +
                                          ti(trend_burnedarea, trend_spei, trend_heatwave, bs=c("ts", "ts", "ts"), k=c(5,5,5)) +
                                          ti(trend_forestpercentage, trend_spei, trend_heatwave, bs=c("ts", "ts", "ts"), k=c(5,5,5)) +
                                          ti(trend_forestpercentage, trend_annualtemp, trend_re, bs=c("ts", "ts", "ts"), k=c(5,5,5)), 
                                        data = forest_df_noNA_nooutliers %>% dplyr::select(-c(ndvi_monthly_lin_dec, 
                                                                                               ndvi_monthly_lin_inc,
                                                                                               ndvi_monthly_step_dec, 
                                                                                               ndvi_monthly_stable,
                                                                                               ndvi_monthly_quad_dec_acc,
                                                                                               ndvi_monthly_quad_dec_decc,
                                                                                               ndvi_monthly_quad_inc_acc,
                                                                                               ndvi_monthly_quad_inc_dec)), 
                                        method = "fREML", 
                                        family = quasibinomial(link = "logit")); Sys.time()
summary(x_model)

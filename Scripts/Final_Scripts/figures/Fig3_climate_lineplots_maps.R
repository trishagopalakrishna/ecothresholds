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
analyses_df_noNA_nooutliers <- read_rds(here("Scripts", "Final_Scripts", "analyses", "analyses_df_noNA_nooutliers.rdata"))

greening_gam <- read_rds(here("Scripts", "Final_Scripts", "analyses", "onemodel_inc.rds"))
browning_gam <- read_rds(here("Scripts", "Final_Scripts", "analyses" ,"onemodel_dec.rds"))

#2. Calculation of univariate smooths and se for greening and browning model
x1 <- -43.2781
y1 <- -4.0213
xaxis_range_len <- 200

grid_creation_prediction <- function (df_analyses, 
                                      xaxis_variable_in_quotes, 
                                      gam_model){
  #Creating new data and running model
  df_select <- df_analyses %>% dplyr::select(-c(notrend))
  
  var_x <- subset( df_select, select = xaxis_variable_in_quotes)
  
  xaxis_range <- seq(min(var_x),max(var_x),len = xaxis_range_len)
  
  df_means<- df_select %>% dplyr::select(-c(x,y,xaxis_variable_in_quotes)) %>%
    summarise(across(where(is.numeric), ~mean(.x, na.rm= T)))
  
  newd <- data.frame(expand.grid(xaxis_range))
  names(newd) <- c(xaxis_variable_in_quotes)
  newd <- newd %>% bind_cols(df_means) %>% mutate(x= x1, y=y1)
  pred <- mgcv::predict.bam(gam_model, newdata = newd, type="terms") 
  pred <- as.data.frame(pred)
  
  #Selecting the coefficients
  x1Eff_inv <- inv.logit(pred[[paste0("s(",xaxis_variable_in_quotes, ")")]] + gam_model$coefficients[1])
  
  #Estimating 95% CI uncertainty
  X <- mgcv::predict.bam(gam_model, newdata = newd, type ="lpmatrix")
  n.sims <- 1000
  b_sims <- rmvn(n.sims,coef(gam_model), gam_model$Vp) #beta unceratinty
  
  smooth_index_trial1 <- grep(paste0("s(",xaxis_variable_in_quotes, ")"),names(coef(gam_model)),fixed = T) #selecting smooths
  smooth_index <- c(smooth_index_trial1)
  
  X1Eff_without_latlong <- tcrossprod( X[,smooth_index] , b_sims[,smooth_index] ) #spline computation
  
  calc_lower_XEff_1<- apply(X1Eff_without_latlong[1:200,], 1, quantile,probs = 0.025)
  invlogit_lower_xEff_1 <- inv.logit(calc_lower_XEff_1 + gam_model$coefficients[1])
  lower_xEff_1<- as.data.frame(invlogit_lower_xEff_1)
  calc_higher_xEff_1 <- apply(X1Eff_without_latlong[1:200,], 1, quantile,probs = 0.975)
  invlogit_higher_xEff_1 <- inv.logit(calc_higher_xEff_1 + gam_model$coefficients[1])
  higher_xEff_1 <- as.data.frame(invlogit_higher_xEff_1)
  ci_xEff <- bind_cols(lower_xEff_1, higher_xEff_1, xaxis_range)
  names(ci_xEff)<- c("lower.ci", "upper.ci", xaxis_variable_in_quotes)
  
  plot_df <- data.frame(xaxis_colname = xaxis_range)
  
  average_plot_df <- plot_df %>% mutate(x1Eff_colname = x1Eff_inv)
  names(average_plot_df)<- c(xaxis_variable_in_quotes, "main_effect")
  average_plot_df <- left_join(average_plot_df, ci_xEff, by = xaxis_variable_in_quotes)
  
  final_plot_df <- bind_rows( average_plot_df)
  final_plot_df 
}

#Mean variables
greening_mean_at <- grid_creation_prediction(analyses_df_noNA_nooutliers, 
                                             "mean_at", greening_gam)
greening_mean_re <- grid_creation_prediction(analyses_df_noNA_nooutliers, 
                                             "mean_re", greening_gam)
greening_mean_spei <- grid_creation_prediction(analyses_df_noNA_nooutliers, 
                                               "mean_spei", greening_gam)
greening_mean_heatwave <- grid_creation_prediction(analyses_df_noNA_nooutliers, 
                                                   "mean_heatwave", greening_gam)

browning_mean_at <- grid_creation_prediction(analyses_df_noNA_nooutliers, 
                                             "mean_at", browning_gam)
browning_mean_re <- grid_creation_prediction(analyses_df_noNA_nooutliers, 
                                             "mean_re", browning_gam)
browning_mean_spei <- grid_creation_prediction(analyses_df_noNA_nooutliers, 
                                               "mean_spei", browning_gam)
browning_mean_heatwave <- grid_creation_prediction(analyses_df_noNA_nooutliers, 
                                                   "mean_heatwave", browning_gam)
#Trend variables
greening_trend_re <- grid_creation_prediction(analyses_df_noNA_nooutliers, 
                                              "trend_re", greening_gam)
greening_trend_at <- grid_creation_prediction(analyses_df_noNA_nooutliers, 
                                              "trend_annualtemp", greening_gam)
greening_trend_spei <- grid_creation_prediction(analyses_df_noNA_nooutliers, 
                                                "trend_spei", greening_gam)
greening_trend_heatwave <- grid_creation_prediction(analyses_df_noNA_nooutliers, 
                                                    "trend_heatwave", greening_gam)

browning_trend_re <- grid_creation_prediction(analyses_df_noNA_nooutliers, 
                                              "trend_re", browning_gam)
browning_trend_at <- grid_creation_prediction(analyses_df_noNA_nooutliers, 
                                              "trend_annualtemp", browning_gam)
browning_trend_spei <- grid_creation_prediction(analyses_df_noNA_nooutliers, 
                                                "trend_spei", browning_gam)
browning_trend_heatwave <- grid_creation_prediction(analyses_df_noNA_nooutliers, 
                                                    "trend_heatwave", browning_gam)

combined_plot_df_prep <- function (greening_results, browning_results){
  greening_results <- greening_results %>% mutate(TrajType = "Greening")
  browning_results <- browning_results %>% mutate(TrajType = "Browning")
  plot_df <- bind_rows (greening_results, browning_results)
  plot_df
}

#Mean variables
mean_at_plot_df <- combined_plot_df_prep(greening_mean_at, browning_mean_at)
mean_re_plot_df <- combined_plot_df_prep(greening_mean_re, browning_mean_re)
mean_spei_plot_df <- combined_plot_df_prep(greening_mean_spei, browning_mean_spei)
mean_heatwave_plot_df <- combined_plot_df_prep(greening_mean_heatwave, browning_mean_heatwave)

#Trend variables
trend_at_plot_df <- combined_plot_df_prep(greening_trend_at, browning_trend_at)
trend_re_plot_df <- combined_plot_df_prep(greening_trend_re, browning_trend_re)
trend_spei_plot_df <- combined_plot_df_prep(greening_trend_spei, browning_trend_spei)
trend_heatwave_plot_df <- combined_plot_df_prep(greening_trend_heatwave, browning_trend_heatwave)

univariate_relationships_list <- list(mean_at_plot_df, mean_re_plot_df, mean_spei_plot_df, mean_heatwave_plot_df,
                                      trend_at_plot_df, trend_re_plot_df, trend_spei_plot_df, trend_heatwave_plot_df)

remove(mean_at_plot_df,mean_re_plot_df, mean_spei_plot_df, mean_heatwave_plot_df, trend_at_plot_df, trend_re_plot_df, trend_spei_plot_df, trend_heatwave_plot_df)

remove(greening_mean_at, greening_mean_re, greening_mean_spei, greening_mean_heatwave, greening_trend_at, greening_trend_re, greening_trend_spei, greening_trend_heatwave, browning_mean_at, browning_mean_re, browning_mean_spei, browning_mean_heatwave, browning_trend_at, browning_trend_re, browning_trend_spei, browning_trend_heatwave)

#3. Line plots of partial effects for only trends climate variables
univariate_plotting <- function (var, x_label, univariate_relationship_df){
  
  data_df <- analyses_df_noNA_nooutliers %>%
    dplyr::select(all_of(var)) %>%
    dplyr::rename(z = 1)
  density <- ggplot(data_df, aes(x = z)) +
    geom_density(colour = "grey70", fill = "grey90", alpha = 0.6) +
    geom_vline(xintercept = quantile(data_df$z, probs=0.1), linetype="dashed") +
    geom_vline(xintercept = quantile(data_df$z, probs=0.9), linetype="dashed") +
    theme_void()
  
  plot <- ggplot(univariate_relationship_df, aes(x = .data[[var]], y = main_effect, group = TrajType, colour = TrajType, fill = TrajType)) +
    geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), alpha = 0.2) + 
    geom_line(alpha = 0.8, lwd = 0.6) +
    coord_cartesian(xlim = c(quantile(univariate_relationship_df$var, probs=0.1),quantile(univariate_relationship_df$var, probs=0.9)), ylim = c(0,1)) +
    geom_hline(yintercept= summary(analyses_df_noNA_nooutliers$increase)[[4]]/25, linetype='longdash', col = 'grey')+
    geom_hline(yintercept= summary(analyses_df_noNA_nooutliers$decrease)[[4]]/25, linetype='longdash', col = 'grey')+
    scale_fill_manual(values = c("#996633", "#669900")) + 
    scale_color_manual(values = c("#663300", "#61A36A")) +
    scale_x_continuous( labels = label_number()) +
    labs(x = x_label,
         y = "Probability of trajectory type") +
    theme_classic(base_size = 16) + theme(legend.position="none")
  
  x<- patchwork::wrap_elements(density) + 
    patchwork::wrap_elements(plot) +
    patchwork::plot_layout(ncol = 1, nrow = 2, widths = 4, heights = c(1, 4))
  x
  
}

trend_at_plot <- univariate_plotting("trend_annualtemp", "trend in annual temperature (degree centigrade/year)", univariate_relationships_list[[5]])
trend_re_plot <- univariate_plotting("trend_re", "trend in relative entropy (per year)", univariate_relationships_list[[6]])
trend_spei_plot <- univariate_plotting("trend_spei", "trend in SPEI (per year)", univariate_relationships_list[[7]])

additional_plots <- function (var, univariate_relationship_df, x_label, fill_color, line_color){
  data_df <- analyses_df_noNA_nooutliers %>%
    dplyr::select(var) %>%
    dplyr::rename(z = 1)
  density <- ggplot(data_df, aes(x = z)) +
    geom_density(colour = "grey70", fill = "grey90", alpha = 0.6) +
    geom_vline(xintercept = quantile(data_df$z, probs=0.1), linetype="dashed") +
    geom_vline(xintercept = quantile(data_df$z, probs=0.9), linetype="dashed") +
    theme_void() 
  
  plot <- ggplot(univariate_relationship_df, aes(x = .data[[var]], y = main_effect, group = TrajType, colour = TrajType, fill = TrajType)) +
    geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), alpha = 0.2) + 
    geom_line(alpha = 0.8, lwd = 0.6) +
    coord_cartesian(xlim = c(quantile(univariate_relationship_df$var0, probs=0.1),quantile(univariate_relationship_df$var0, probs=0.9)), ylim = c(0,1)) +
    geom_hline(yintercept= summary(analyses_df_noNA_nooutliers$increase)[[4]]/25, linetype='longdash', col = 'grey')+
    scale_fill_manual(values = fill_color) + 
    scale_color_manual(values = line_color) +
    scale_x_continuous( labels = label_number()) +
    labs(x = x_label,
         y = "Probability of trajectory type") +
    theme_classic(base_size = 16) + theme(legend.position="none")
  
  x <- density + plot +
    patchwork::plot_layout(ncol = 1, nrow = 2, widths = 4, heights = c(1, 4))
  x
}
trend_heatwave_greening <- additional_plots("trend_heatwave", 
                                            univariate_relationships_list[[8]] %>% filter(TrajType=="Greening"), 
                                            "trend in heatwaves","#669900", "#61A36A")


#4. Greening and browning probability maps of above trend variables
#4 (i) Estimating map of greening and browning of above trend climate variables
mapping_probability_function <- function (driver_variable, gam_model){
  df_means <- analyses_df_noNA_nooutliers %>% 
    dplyr::select(-c(x, y, driver_variable)) %>% 
    summarise(across(where(is.numeric), ~mean(.x, na.rm= T)))
  
  df_to_predict<- analyses_df_noNA_nooutliers %>% 
    dplyr::select(c(x, y, driver_variable)) %>% 
    bind_cols(df_means)
  
  pred <- mgcv::predict.bam(gam_model, newdata = df_to_predict, type="terms") 
  pred <- as.data.frame(pred)
  pred <- pred %>% bind_cols(analyses_df_noNA_nooutliers %>% 
                               dplyr::select(c(x, y, driver_variable)))
  
  coef_index_trial1 <- grep(paste0("s(",driver_variable,")"),names(pred),fixed = T)
  coef_index_trial2 <- grep(paste0("ti(mean_at,x,y)"),names(pred),fixed = T)
  coef_index_trial3 <- grep(paste0("ti(trend_annualtemp,x,y)"),names(pred),fixed = T)
  
  effects_cols_indices <- c(coef_index_trial1, coef_index_trial2, coef_index_trial3)
  
  x1Eff <- rowSums(pred [effects_cols_indices])
  x1Eff_inv <- inv.logit(x1Eff + gam_model$coefficients[1])
  x1Eff_inv <- as.data.frame(x1Eff_inv)
  x1Eff_inv <- x1Eff_inv %>% bind_cols(analyses_df_noNA_nooutliers %>% 
                                         dplyr::select(c(x, y)))
  names(x1Eff_inv)[1]<- "main_effect"
  x1Eff_inv <- x1Eff_inv %>% dplyr::select(x,y, main_effect)
  x_rast <- rast(as.data.frame(x1Eff_inv %>% drop_na()), type="xyz")
  
  x_rast
}

#Greening
Sys.time();greening_trend_at <- mapping_probability_function("trend_annualtemp", greening_gam); Sys.time()
Sys.time();greening_trend_re <-  mapping_probability_function("trend_re", greening_gam); Sys.time()
Sys.time();greening_trend_heatwave <- mapping_probability_function("trend_heatwave", greening_gam); Sys.time()
Sys.time();greening_trend_spei <-  mapping_probability_function("trend_spei", greening_gam); Sys.time()

#Browning
Sys.time();browning_trend_at <- mapping_probability_function("trend_annualtemp", browning_gam); Sys.time()
Sys.time();browning_trend_re <-  mapping_probability_function("trend_re", browning_gam); Sys.time()
Sys.time();browning_trend_spei <-  mapping_probability_function("trend_spei", browning_gam); Sys.time()

#4 (ii) Maps of greening and browning of above trend climate variables
d_trans<- st_read(here ("Data", "Cattelanetal_Clustering", "cerrado_climate_zones.shp"))
d_trans<- st_transform(d_trans, crs = 4326)
d_trans<- d_trans %>% mutate(area_km= terra::expanse(vect(d_trans), unit="km"))
cerrado<- d_trans %>% st_union()

aggregate_mapping_function <- function (partialeffect_raster, greens_or_brown_values){
  x_map<- tm_shape(cerrado)+ tm_borders()+
    tm_polygons(fill.scale = tm_scale_categorical(n=1, values = "grey"))+
    tm_shape (partialeffect_raster) + 
    tm_raster(col.scale = tm_scale_continuous(values = greens_or_brown_values),
              col.legend = tm_legend(show= FALSE, title = "Probability")) + 
    tm_layout(frame = FALSE)
  x_map
}

Sys.time(); map_green_trend_at <- aggregate_mapping_function(greening_trend_at, "greens"); Sys.time()
Sys.time(); map_green_trend_re <- aggregate_mapping_function(greening_trend_re, "greens"); Sys.time()
Sys.time(); map_green_trend_heatwave <- aggregate_mapping_function(greening_trend_heatwave, "greens"); Sys.time()
Sys.time(); map_green_trend_spei <- aggregate_mapping_function(greening_trend_spei, "greens"); Sys.time()

Sys.time(); map_brown_trend_at <- aggregate_mapping_function(browning_trend_at, "brown"); Sys.time()
Sys.time(); map_brown_trend_re <- aggregate_mapping_function(browning_trend_re, "brown"); Sys.time()
Sys.time(); map_brown_trend_spei <- aggregate_mapping_function(browning_trend_spei, "brown"); Sys.time()

green_legend_only <- tm_shape(cerrado)+ tm_borders()+
  tm_polygons(fill.scale = tm_scale_categorical(n=1, values = "grey"))+
  tm_shape (greening_trend_at) + 
  tm_raster(col.scale = tm_scale_continuous(values = "greens"),
            col.legend = tm_legend(title = "Probability", orientation = "landscape",
                                   reverse = T,
                                   frame = FALSE)) + 
  tm_layout( legend.width = 20,
    legend.only = TRUE, frame = FALSE)
tmap_save(green_legend_only, here("Outputs", "TrendsResults", "Figures", "green_legend.png"),
          height = 5, width = 8, units = "cm", dpi=700)

brown_legend_only <- tm_shape(cerrado)+ tm_borders()+
  tm_polygons(fill.scale = tm_scale_categorical(n=1, values = "grey"))+
  tm_shape (browning_trend_at) + 
  tm_raster(col.scale = tm_scale_continuous(values = "brown"),
            col.legend = tm_legend(title = "Probability", orientation = "landscape",
                                   reverse = T,
                                   frame = FALSE)) + 
  tm_layout(legend.width = 20, legend.only = TRUE, frame = FALSE)
tmap_save(brown_legend_only, here("Outputs", "TrendsResults", "Figures", "brown_legend.png"),
          height = 5, width = 8, units = "cm", dpi=700)

#5 Combining plots
greening_maps <- tmap_arrange(map_green_trend_at, map_green_trend_re, 
                              map_green_trend_spei, map_green_trend_heatwave,
                              ncol = 1, nrow = 4)

browning_maps <- tmap_arrange(map_brown_trend_at, map_brown_trend_re, 
                              map_brown_trend_spei, 
                              ncol = 1, nrow = 3)

line_plots <- ggpubr::ggarrange(plotlist = list(trend_at_plot, trend_re_plot, 
                                                trend_spei_plot,
                                               trend_heatwave_greening),
                               labels=c("a)", "b)", "c)", "d)"),
                               ncol = 1, nrow = 4)


library(patchwork)

at_plot <- patchwork::wrap_elements(trend_at_plot) + 
  patchwork::wrap_elements(tmap_grob(map_green_trend_at)) +
  patchwork::wrap_elements(tmap_grob(map_brown_trend_at)) +
  patchwork::plot_layout(ncol = 3, nrow = 1) +
  patchwork::plot_annotation(title = "a)")

re_plot <-  patchwork::wrap_elements(trend_re_plot) +
  patchwork::wrap_elements(tmap_grob(map_green_trend_re)) +
  patchwork::wrap_elements(tmap_grob(map_brown_trend_re)) +
  patchwork::plot_layout(ncol = 3, nrow = 1)+
  patchwork::plot_annotation(title = "b)")

spei_plot <-  patchwork::wrap_elements(trend_spei_plot) +
  patchwork::wrap_elements(tmap_grob(map_green_trend_spei)) +
  patchwork::wrap_elements(tmap_grob(map_brown_trend_spei)) +
  patchwork::plot_layout(ncol = 3, nrow = 1)+
  patchwork::plot_annotation(title = "c)")
  
heatwave_plot<-  patchwork::wrap_elements(trend_heatwave_greening) +
  patchwork::wrap_elements(tmap_grob(map_green_trend_heatwave)) +
    patchwork::plot_layout(ncol = 3, nrow = 1)+
  patchwork::plot_annotation(title = "d)")

trial <- at_plot / re_plot / spei_plot / heatwave_plot
ggsave(here("Outputs", "TrendsResults", "Figures", "trial.png"),
       trial, dpi = 700, height = 55, width= 55, units = "cm")


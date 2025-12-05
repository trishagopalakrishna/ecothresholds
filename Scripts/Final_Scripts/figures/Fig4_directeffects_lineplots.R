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
analyses_df_noNA_nooutliers <- read_rds(here("Scripts", "Final_Scripts", "analyses", "analyses_df_noNA_nooutliers.rds"))

greening_gam <- read_rds(here("Scripts", "Final_Scripts", "analyses", "onemodel_inc.rds"))

#2. Calculation of smooths and se for greening and browning model only for direct causal relationships
x1 <- -43.2781
y1 <- -4.0213
xaxis_range_len <- 200

#2. (i) For composite climate variables ti(trend_at, trend_re) 
composite_grid_creation_prediction1 <- function (df_analyses, 
                                            xaxis_variable_in_quotes, 
                                            line_variable_in_quotes, 
                                            gam_model){
  #Creating new data and running model
  df_select <- df_analyses %>% dplyr::select(-c(notrend))
  
  var_x <- subset( df_select, select = xaxis_variable_in_quotes)
  var_line <- subset( df_select, select = line_variable_in_quotes)
  var_line <- as.vector(var_line) %>% unlist()
  
  xaxis_range <- seq(min(var_x),max(var_x),len = xaxis_range_len)
  line_range <- c((mean(var_line) - sd(var_line)), 0 , (mean(var_line) + sd(var_line))) #changed middle scenario to be 0, so no change in trend
  
  df_means<- df_select %>% dplyr::select(-c(x,y,xaxis_variable_in_quotes, line_variable_in_quotes)) %>%
    summarise(across(where(is.numeric), ~mean(.x, na.rm= T)))
  
  newd <- data.frame(expand.grid(xaxis_range, line_range))
  names(newd) <- c(xaxis_variable_in_quotes, line_variable_in_quotes)
  newd <- newd %>% bind_cols(df_means) %>% mutate(x= x1, y=y1)
  pred <- mgcv::predict.bam(gam_model, newdata = newd, type="terms") 
  pred <- as.data.frame(pred)
  
  #Selecting the coefficients
  ind <- 1:xaxis_range_len
  ind2 <- (xaxis_range_len + 1):(xaxis_range_len*2)
  ind3 <- ((xaxis_range_len*2)+1):(xaxis_range_len*3)
  
  xy_interaction_effects <- grep("x,y",names(pred),fixed = T)
  pred_without_xy_interaction_effects <- (subset(pred, select = -xy_interaction_effects))
  
  coef_index_trial1 <- grep(paste0("s(",xaxis_variable_in_quotes, ")"),names(pred_without_xy_interaction_effects),fixed = T )
  coef_index_trial2 <- grep(paste0("s(",line_variable_in_quotes, ")"),names(pred_without_xy_interaction_effects),fixed = T )
  coef_index_trial3 <- grep(paste0("ti(", xaxis_variable_in_quotes, ",", line_variable_in_quotes, ")"),names(pred_without_xy_interaction_effects),fixed = T)
  
  effects_cols_indices <- c(coef_index_trial1, coef_index_trial2, coef_index_trial3)
  
  x1Eff <- rowSums(pred_without_xy_interaction_effects [ind, effects_cols_indices])
  x1Eff_inv <- inv.logit(x1Eff + gam_model$coefficients[1])
  x2Eff <- rowSums(pred_without_xy_interaction_effects [ind2, effects_cols_indices])
  x2Eff_inv <- inv.logit(x2Eff + gam_model$coefficients[1])
  x3Eff <- rowSums(pred_without_xy_interaction_effects [ind3, effects_cols_indices])
  x3Eff_inv <- inv.logit(x3Eff + gam_model$coefficients[1])
  
  #Estimating 95% CI uncertainty
  X <- mgcv::predict.bam(gam_model, newdata = newd, type ="lpmatrix")
  n.sims <- 1000
  b_sims <- rmvn(n.sims,coef(gam_model), gam_model$Vp) #beta unceratinty
  
  smooth_index_trial1 <- grep(paste0("s(",xaxis_variable_in_quotes, ")"),names(coef(gam_model)),fixed = T) #selecting smooths
  smooth_index_trial2 <- grep(paste0("s(",line_variable_in_quotes, ")"),names(coef(gam_model)),fixed = T )
  smooth_index_trial3 <- grep(paste0("ti(", xaxis_variable_in_quotes, ",", line_variable_in_quotes, ")"),names(coef(gam_model)),fixed = T )
  
  smooth_index <- c(smooth_index_trial1, smooth_index_trial2, smooth_index_trial3)
  
  X1Eff_without_latlong <- tcrossprod( X[,smooth_index] , b_sims[,smooth_index] ) #spline computation
  
  calc_lower_XEff_1<- apply(X1Eff_without_latlong[1:200,], 1, quantile,probs = 0.025)
  invlogit_lower_xEff_1 <- inv.logit(calc_lower_XEff_1 + gam_model$coefficients[1])
  lower_xEff_1<- as.data.frame(invlogit_lower_xEff_1)
  calc_higher_xEff_1 <- apply(X1Eff_without_latlong[1:200,], 1, quantile,probs = 0.975)
  invlogit_higher_xEff_1 <- inv.logit(calc_higher_xEff_1 + gam_model$coefficients[1])
  higher_xEff_1 <- as.data.frame(invlogit_higher_xEff_1)
  ci_xEff <- bind_cols(lower_xEff_1, higher_xEff_1, xaxis_range)
  names(ci_xEff)<- c("lower.ci", "upper.ci", xaxis_variable_in_quotes)
  
  calc_lower_XEff_2<- apply(X1Eff_without_latlong[201:400,], 1, quantile,probs = 0.025)
  invlogit_lower_xEff_2 <- inv.logit(calc_lower_XEff_2 + gam_model$coefficients[1])
  lower_xEff_2 <- as.data.frame(invlogit_lower_xEff_2)
  calc_higher_xEff_2 <- apply(X1Eff_without_latlong[201:400,], 1, quantile,probs = 0.975)
  invlogit_higher_xEff_2 <- inv.logit(calc_higher_xEff_2 + gam_model$coefficients[1])
  higher_xEff_2 <- as.data.frame(invlogit_higher_xEff_2)
  ci_xEff_2 <- bind_cols(lower_xEff_2, higher_xEff_2, xaxis_range)
  names(ci_xEff_2)<- c("lower.ci", "upper.ci", xaxis_variable_in_quotes)
  
  calc_lower_XEff_3 <- apply(X1Eff_without_latlong[401:600,], 1, quantile,probs = 0.025)
  invlogit_lower_xEff_3 <- inv.logit(calc_lower_XEff_3 + gam_model$coefficients[1])
  lower_xEff_3 <- as.data.frame(invlogit_lower_xEff_3)
  calc_higher_xEff_3 <- apply(X1Eff_without_latlong[401:600,], 1, quantile,probs = 0.975)
  invlogit_higher_xEff_3 <- inv.logit(calc_higher_xEff_3 + gam_model$coefficients[1])
  higher_xEff_3 <- as.data.frame(invlogit_higher_xEff_3)
  ci_xEff_3 <- bind_cols(lower_xEff_3, higher_xEff_3, xaxis_range)
  names(ci_xEff_3)<- c("lower.ci", "upper.ci", xaxis_variable_in_quotes)
  
  plot_df <- data.frame(xaxis_colname = xaxis_range)
  low_plot_df <- plot_df %>% mutate(x1Eff_colname = x1Eff_inv, line_variable_type = line_range[[1]])
  names(low_plot_df)<- c(xaxis_variable_in_quotes, "main_effect", "line_variable_type")
  low_plot_df <- left_join(low_plot_df, ci_xEff, by = xaxis_variable_in_quotes)
  low_plot_df <- low_plot_df %>% mutate(line_variable_type2 = paste0("decreasing (",round(line_variable_type, digits =4), ")"))
  
  average_plot_df <- plot_df %>% mutate(x1Eff_colname = x2Eff_inv, line_variable_type = line_range[[2]])
  names(average_plot_df)<- c(xaxis_variable_in_quotes, "main_effect", "line_variable_type")
  average_plot_df <- left_join(average_plot_df, ci_xEff_2, by = xaxis_variable_in_quotes)
  average_plot_df <- average_plot_df %>% mutate(line_variable_type2 = paste0("no change (",round(line_variable_type, digits =4), ")"))
  
  high_plot_df <- plot_df %>% mutate(x1Eff_colname = x3Eff_inv, line_variable_type = line_range[[3]])
  names(high_plot_df)<- c(xaxis_variable_in_quotes, "main_effect", "line_variable_type")
  high_plot_df <- left_join(high_plot_df, ci_xEff_3, by = xaxis_variable_in_quotes)
  high_plot_df <- high_plot_df %>% mutate(line_variable_type2 = paste0("increasing (",round(line_variable_type, digits =4), ")"))
  
  final_plot_df <- bind_rows(low_plot_df, average_plot_df, high_plot_df)
  final_plot_df 
}

greening_trend_at_re <- composite_grid_creation_prediction1(analyses_df_noNA_nooutliers, 
                                                       "trend_annualtemp", 
                                                       "trend_re",
                                                       greening_gam)

# & ti(trend_heatwaves, trend_spei). Considering only one "more drought" condition value
composite_grid_creation_prediction2 <- function (df_analyses, 
                                                 xaxis_variable_in_quotes, 
                                                 line_variable_in_quotes, 
                                                 gam_model){
  #Creating new data and running model
  df_select <- df_analyses %>% dplyr::select(-c(notrend))
  
  var_x <- subset( df_select, select = xaxis_variable_in_quotes)
  var_line <- subset( df_select, select = line_variable_in_quotes)
  var_line <- as.vector(var_line) %>% unlist()
  
  xaxis_range <- seq(min(var_x),max(var_x),len = xaxis_range_len)
  line_range <- c((mean(var_line) - sd(var_line)), 0) #changed middle scenario to be 0, so no change in trend
  
  df_means<- df_select %>% dplyr::select(-c(x,y,xaxis_variable_in_quotes, line_variable_in_quotes)) %>%
    summarise(across(where(is.numeric), ~mean(.x, na.rm= T)))
  
  newd <- data.frame(expand.grid(xaxis_range, line_range))
  names(newd) <- c(xaxis_variable_in_quotes, line_variable_in_quotes)
  newd <- newd %>% bind_cols(df_means) %>% mutate(x= x1, y=y1)
  pred <- mgcv::predict.bam(gam_model, newdata = newd, type="terms") 
  pred <- as.data.frame(pred)
  
  #Selecting the coefficients
  ind <- 1:xaxis_range_len
  ind2 <- (xaxis_range_len + 1):(xaxis_range_len*2)
  
  xy_interaction_effects <- grep("x,y",names(pred),fixed = T)
  pred_without_xy_interaction_effects <- (subset(pred, select = -xy_interaction_effects))
  
  coef_index_trial1 <- grep(paste0("s(",xaxis_variable_in_quotes, ")"),names(pred_without_xy_interaction_effects),fixed = T )
  coef_index_trial2 <- grep(paste0("s(",line_variable_in_quotes, ")"),names(pred_without_xy_interaction_effects),fixed = T )
  coef_index_trial3 <- grep(paste0("ti(", xaxis_variable_in_quotes, ",", line_variable_in_quotes, ")"),names(pred_without_xy_interaction_effects),fixed = T)
  
  effects_cols_indices <- c(coef_index_trial1, coef_index_trial2, coef_index_trial3)
  
  x1Eff <- rowSums(pred_without_xy_interaction_effects [ind, effects_cols_indices])
  x1Eff_inv <- inv.logit(x1Eff + gam_model$coefficients[1])
  x2Eff <- rowSums(pred_without_xy_interaction_effects [ind2, effects_cols_indices])
  x2Eff_inv <- inv.logit(x2Eff + gam_model$coefficients[1])
  
  #Estimating 95% CI uncertainty
  X <- mgcv::predict.bam(gam_model, newdata = newd, type ="lpmatrix")
  n.sims <- 1000
  b_sims <- rmvn(n.sims,coef(gam_model), gam_model$Vp) #beta unceratinty
  
  smooth_index_trial1 <- grep(paste0("s(",xaxis_variable_in_quotes, ")"),names(coef(gam_model)),fixed = T) #selecting smooths
  smooth_index_trial2 <- grep(paste0("s(",line_variable_in_quotes, ")"),names(coef(gam_model)),fixed = T )
  smooth_index_trial3 <- grep(paste0("ti(", xaxis_variable_in_quotes, ",", line_variable_in_quotes, ")"),names(coef(gam_model)),fixed = T )
  
  smooth_index <- c(smooth_index_trial1, smooth_index_trial2, smooth_index_trial3)
  
  X1Eff_without_latlong <- tcrossprod( X[,smooth_index] , b_sims[,smooth_index] ) #spline computation
  
  calc_lower_XEff_1<- apply(X1Eff_without_latlong[1:200,], 1, quantile,probs = 0.025)
  invlogit_lower_xEff_1 <- inv.logit(calc_lower_XEff_1 + gam_model$coefficients[1])
  lower_xEff_1<- as.data.frame(invlogit_lower_xEff_1)
  calc_higher_xEff_1 <- apply(X1Eff_without_latlong[1:200,], 1, quantile,probs = 0.975)
  invlogit_higher_xEff_1 <- inv.logit(calc_higher_xEff_1 + gam_model$coefficients[1])
  higher_xEff_1 <- as.data.frame(invlogit_higher_xEff_1)
  ci_xEff <- bind_cols(lower_xEff_1, higher_xEff_1, xaxis_range)
  names(ci_xEff)<- c("lower.ci", "upper.ci", xaxis_variable_in_quotes)
  
  calc_lower_XEff_2<- apply(X1Eff_without_latlong[201:400,], 1, quantile,probs = 0.025)
  invlogit_lower_xEff_2 <- inv.logit(calc_lower_XEff_2 + gam_model$coefficients[1])
  lower_xEff_2 <- as.data.frame(invlogit_lower_xEff_2)
  calc_higher_xEff_2 <- apply(X1Eff_without_latlong[201:400,], 1, quantile,probs = 0.975)
  invlogit_higher_xEff_2 <- inv.logit(calc_higher_xEff_2 + gam_model$coefficients[1])
  higher_xEff_2 <- as.data.frame(invlogit_higher_xEff_2)
  ci_xEff_2 <- bind_cols(lower_xEff_2, higher_xEff_2, xaxis_range)
  names(ci_xEff_2)<- c("lower.ci", "upper.ci", xaxis_variable_in_quotes)

  
  plot_df <- data.frame(xaxis_colname = xaxis_range)
  low_plot_df <- plot_df %>% mutate(x1Eff_colname = x1Eff_inv, line_variable_type = line_range[[1]])
  names(low_plot_df)<- c(xaxis_variable_in_quotes, "main_effect", "line_variable_type")
  low_plot_df <- left_join(low_plot_df, ci_xEff, by = xaxis_variable_in_quotes)
  low_plot_df <- low_plot_df %>% mutate(line_variable_type2 = paste0("decreasing (",round(line_variable_type, digits =4), ")"))
  
  average_plot_df <- plot_df %>% mutate(x1Eff_colname = x2Eff_inv, line_variable_type = line_range[[2]])
  names(average_plot_df)<- c(xaxis_variable_in_quotes, "main_effect", "line_variable_type")
  average_plot_df <- left_join(average_plot_df, ci_xEff_2, by = xaxis_variable_in_quotes)
  average_plot_df <- average_plot_df %>% mutate(line_variable_type2 = paste0("no change (",round(line_variable_type, digits =4), ")"))
  
  final_plot_df <- bind_rows(low_plot_df, average_plot_df)
  final_plot_df 
}

greening_trend_heatwave_spei <- composite_grid_creation_prediction2(analyses_df_noNA_nooutliers, 
                                                               "trend_heatwave", 
                                                               "trend_spei",
                                                               greening_gam)

composite_combined_plot_df_prep <- function (greening_results, browning_results){
  greening_results <- greening_results %>% mutate(TrajType = "Greening")
  greening_results
}

composite_at_re_plot_df <- composite_combined_plot_df_prep(greening_trend_at_re)
composite_heatwave_spei_plot_df <- composite_combined_plot_df_prep(greening_trend_heatwave_spei)

composite_plot_df_list <- list(composite_at_re_plot_df, composite_heatwave_spei_plot_df)
remove(greening_trend_at_re, greening_trend_heatwave_spei)
remove(composite_combined_plot_df_prep, composite_grid_creation_prediction1, composite_grid_creation_prediction2)
remove(composite_at_re_plot_df, composite_heatwave_spei_plot_df)

#2. (ii) For univariate smooth relationships of anthropic pressures, fire activity and savanna vegetation cover

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

greening_anthropicdist <- grid_creation_prediction(analyses_df_noNA_nooutliers, 
                                                   "anthropic_dist", greening_gam)
greening_trend_burnedarea <- grid_creation_prediction(analyses_df_noNA_nooutliers, 
                                                     "trend_burnedarea", greening_gam)
greening_mean_savanna <- grid_creation_prediction(analyses_df_noNA_nooutliers, 
                                                 "mean_savanna_percentage", greening_gam)

combined_plot_df_prep <- function (greening_results){
  greening_results <- greening_results %>% mutate(TrajType = "Greening")
  greening_results
}

anthropic_plot_df <- combined_plot_df_prep(greening_anthropicdist)
burnedarea_plot_df <- combined_plot_df_prep(greening_trend_burnedarea)
greening_mean_savanna <- combined_plot_df_prep(greening_mean_savanna) 

univariate_plot_df_list <- list(anthropic_plot_df, burnedarea_plot_df, greening_mean_savanna)

remove(greening_anthropicdist, greening_trend_burnedarea, greening_mean_savanna)
remove(grid_creation_prediction, combined_plot_df_prep)
remove(anthropic_plot_df, burnedarea_plot_df)

#3. Line plots of partial effects/probability for direct causal relationships
greening_intercept_transposed <- inv.logit(greening_gam$coefficients[[1]])

composite_plot_function <- function (interaction_trend_df, xaxis_variable, xaxis_label){
  x_plot<- ggplot(interaction_trend_df, aes(x = .data[[xaxis_variable]], y = main_effect,
                                            color = TrajType, fill =TrajType,
                                            group = interaction(line_variable_type2, TrajType))) +
    geom_line(aes(color = TrajType), linewidth = 1) + 
    geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), alpha = 0.1) +
    coord_cartesian(xlim = c(quantile(interaction_trend_df[[xaxis_variable]], probs=0.1),quantile(interaction_trend_df[[xaxis_variable]], probs=0.9)), ylim = c(0,1))+
    geom_hline(yintercept= greening_intercept_transposed, linetype='longdash', col = 'grey')+
    facet_wrap(~factor(line_variable_type2, unique(interaction_trend_df$line_variable_type2)), nrow=1)+
    theme_classic(base_size = 10) +
    scale_fill_manual(values = c( "#669900"))+
    scale_color_manual(values=c( "#61A36A")) + 
    ylab ("Probability of trajectory type") + 
    xlab(xaxis_label)+
    theme(legend.position="none")
  x_plot
}
composite_plot_trend_at_re <- composite_plot_function(composite_plot_df_list[[1]], "trend_annualtemp", 
                                                      "trend in annual temperature (degree centigrade/year)")
ggsave(here("Outputs", "TrendsResults", "Figures", "Fig3", "panel_re_at2.png"),
       composite_plot_trend_at_re,
       dpi=700, width =15, height = 5, units = "cm")

composite_plot_trend_heatwave_spei <- composite_plot_function(composite_plot_df_list[[2]], "trend_heatwave", 
                                                              "trend in number of heatwaves (per year)")
ggsave(here("Outputs", "TrendsResults", "Figures", "Fig3", "panel_heatwave_spei2.png"),
       composite_plot_trend_heatwave_spei,
       dpi=700, width =15, height = 5, units = "cm")

x_axis_variable_density_plot <- function (xaxis_variable, xaxis_label){
  data_df <- analyses_df_noNA_nooutliers %>%
    dplyr::select(all_of(xaxis_variable)) %>%
    dplyr::rename(z = 1) %>%
    filter(if_all(where(is.numeric),
                  ~ between(.x, quantile(.x, .1), quantile(.x, .9))))
  
  density <- ggplot(data_df, aes(x = z)) +
    geom_density(colour = "grey70", fill = "grey90", alpha = 0.6) +
    theme_classic(base_size = 10)+
    ylab("kernel density") + xlab(xaxis_label)
}
at_density <- x_axis_variable_density_plot("trend_annualtemp", "trend in annual temperature (degree centigrade/year)")
ggsave(here("Outputs", "TrendsResults", "Figures", "Fig3", "at_distribution2.png"),
       at_density,
       dpi=700, width =15, height = 3, units = "cm")
heatwave_density <- x_axis_variable_density_plot("trend_heatwave", "trend in number of heatwaves (per year)")
ggsave(here("Outputs", "TrendsResults", "Figures", "Fig3", "heatwave_distribution2.png"),
       heatwave_density,
       dpi=700, width =15, height = 3, units = "cm")

panel_density_plot1 <- function (panel_term, panel_term_label){
  data_df2 <- analyses_df_noNA_nooutliers %>%
    dplyr::select(all_of(panel_term)) %>%
    dplyr::rename(z = 1)
  
  increasing_value <-  mean(data_df2$z) + sd(data_df2$z)
  decreasing_value <-  mean(data_df2$z) - sd(data_df2$z)
  
  data_df2 <- analyses_df_noNA_nooutliers %>%
  dplyr::select(all_of(panel_term)) %>%
  dplyr::rename(z = 1) %>%
    filter(if_all(where(is.numeric),
                  ~ between(.x, quantile(.x, .1), quantile(.x, .9))))
  
  density2 <- ggplot(data_df2, aes(x = z)) +
    geom_density(colour = "grey70", fill = "grey90", alpha = 0.6) +
    geom_vline(xintercept = increasing_value, linetype="dashed") +
    geom_vline(xintercept = 0, linetype="dashed") +
    geom_vline(xintercept = decreasing_value, linetype="dashed") +
    theme_classic(base_size = 10) +
    ylab("kernel density") + xlab(panel_term_label)+ coord_flip()
}
re_density <- panel_density_plot1("trend_re", "trend in relative entropy \n (per year)")
ggsave(here("Outputs", "TrendsResults", "Figures", "Fig3", "re_distribution2.png"),
       re_density,
       dpi=700, width =5, height = 5, units = "cm")

panel_density_plot2 <- function (panel_term, panel_term_label){
  data_df2 <- analyses_df_noNA_nooutliers %>%
    dplyr::select(all_of(panel_term)) %>%
    dplyr::rename(z = 1)
  
  decreasing_value <-  mean(data_df2$z) - sd(data_df2$z)
  
  data_df2 <- analyses_df_noNA_nooutliers %>%
    dplyr::select(all_of(panel_term)) %>%
    dplyr::rename(z = 1) %>%
    filter(if_all(where(is.numeric),
                  ~ between(.x, quantile(.x, .1), quantile(.x, .9))))
  
  density2 <- ggplot(data_df2, aes(x = z)) +
    geom_density(colour = "grey70", fill = "grey90", alpha = 0.6) +
    geom_vline(xintercept = 0, linetype="dashed") +
    geom_vline(xintercept = decreasing_value, linetype="dashed") +
    theme_classic(base_size = 10) +
    ylab("kernel density") + xlab(panel_term_label)+ coord_flip()
}


spei_density <- panel_density_plot2 ("trend_spei", "trend in spei \n (per year)")
ggsave(here("Outputs", "TrendsResults", "Figures", "Fig3", "spei_distribution2.png"),
       spei_density,
       dpi=700, width = 5, height = 5, units = "cm")

univariate_plotting <- function (var, x_label, univariate_relationship_df){
  
  data_df <- analyses_df_noNA_nooutliers %>%
    dplyr::select(all_of(var)) %>%
    dplyr::rename(z = 1) %>%
    filter(if_all(where(is.numeric),
                  ~ between(.x, quantile(.x, .1), quantile(.x, .9))))
  
  density <- ggplot(data_df, aes(x = z)) +
    geom_density(colour = "grey70", fill = "grey90", alpha = 0.6) +
    theme_void(base_size = 12)
  
  plot <- ggplot(univariate_relationship_df, aes(x = .data[[var]], y = main_effect, group = TrajType, colour = TrajType, fill = TrajType)) +
    geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), alpha = 0.2) + 
    geom_line(alpha = 0.8, lwd = 0.6) +
    coord_cartesian(xlim = c(quantile(univariate_relationship_df$var, probs=0.1),quantile(univariate_relationship_df$var, probs=0.9)), ylim = c(0,1)) +
    geom_hline(yintercept= greening_intercept_transposed, linetype='longdash', col = 'grey')+
    scale_fill_manual(values = c("#669900")) + 
    scale_color_manual(values = c( "#61A36A")) +
    scale_x_continuous( labels = label_number()) +
    labs(x = x_label,
         y = "Probability of trajectory type") +
    theme_classic(base_size = 12) + theme(legend.position="none")
  
  x<- patchwork::wrap_elements(density) + 
    patchwork::wrap_elements(plot) +
    patchwork::plot_layout(ncol = 1, nrow = 2, widths = 4, heights = c(1, 4))
  x
  
}

anthropic_plot <- univariate_plotting("anthropic_dist", "distance to anthropic land use (m)", univariate_plot_df_list[[1]])
ggsave(here("Outputs", "TrendsResults", "Figures", "Fig3", "univariate_anthropicdist2.png"),
       anthropic_plot,
       dpi=700, width =10, height = 10, units = "cm")
burnedarea_plot <- univariate_plotting("trend_burnedarea", "trend in burned area (% per year)", univariate_plot_df_list[[2]])
ggsave(here("Outputs", "TrendsResults", "Figures", "Fig3", "univariate_burnedarea2.png"),
       burnedarea_plot,
       dpi=700, width =10, height = 10, units = "cm")
savanna_plot <- univariate_plotting("mean_savanna_percentage", "average proportion of savanna area (%)", univariate_plot_df_list[[3]])
ggsave(here("Outputs", "TrendsResults", "Figures", "Fig3", "univariate_savanna2.png"),
       savanna_plot,
       dpi=700, width =10, height = 10, units = "cm")


#composite legend only
green_legend_map<- tm_shape(cerrado)+ tm_borders()+
  tm_polygons(fill.scale = tm_scale_categorical(n=1, values = "grey"))+
  tm_shape (greening_trend_at_re) + 
  tm_raster(col.scale = tm_scale_continuous(values = "greens"),
            col.legend = tm_legend(show= TRUE, title = "Probability", reverse = T, frame.color = NA) )+
  tm_layout(legend.only = T)
tmap_save(green_legend_map, here("Outputs", "TrendsResults", "Figures", "Fig3", "green_legend.png"), 
          width = 5, height = 5, units = "cm",
          dpi = 300)



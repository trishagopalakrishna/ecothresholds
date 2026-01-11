require(mgcv)
library(dplyr)
library(tidyverse)
library(tidyr)

browning_gam <- read_rds("onemodel_dec.rds") #load browning gam model with your file path
greening_gam <- read_rds("onemodel_inc.rds") #load greening gam model with your file path
analyses_df <-read_r("analyses_df_noNA_nooutliers.rds") #load data 

summary(analyses_df)
#Note that the "increase" and "decrease" coulmns are the count of the number
#of 1x1km pixels within a 5x5km pixel that has greened or browned. Hence, the 
#maximum value can be 25, which would mean that the entire 5x5km pixel has greened
#or browned. 

#See the mean values in the "increase" and "decrease" coloumns. In the "increase" column,
#18.87/25 ~75% and in the "decrease" column 5.07/25 ~20% which is the average area of a 
#5 x 5km pixel that has greened or browned, considering all 5 x 5km pixels. Hence, the
#75% and 20% are the "reference" or "baseline" proportions of a 5 x 5km pixel that I have tried
#to depict in the plots using dashed grey horizontal lines

#Example line plot code below considering x axis variable to be trend in burned area

#Creating new data grid to predict the model and predicting the browning model
x1 <- -43.2781 #example longitude
y1 <- -4.0213 #example latitude
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
greening_trend_burnedarea <- grid_creation_prediction(analyses_df, 
                                                      "trend_burnedarea", greening_gam)
browning_trend_burnedarea <- grid_creation_prediction(analyses_df, 
                                                      "trend_burnedarea", browning_gam)

combined_plot_df_prep <- function (greening_results, browning_results){
  greening_results <- greening_results %>% mutate(TrajType = "Greening")
  browning_results <- browning_results %>% mutate(TrajType = "Browning")
  plot_df <- bind_rows (greening_results, browning_results)
  plot_df
}
burnedarea_plot_df <- combined_plot_df_prep(greening_trend_burnedarea, browning_trend_burnedarea)

plotting_function <- function (var, x_label, plot_df){
  
  data_df <- analyses_df %>%
    dplyr::select(all_of(var)) %>%
    dplyr::rename(z = 1)
  density <- ggplot(data_df, aes(x = z)) +
    geom_density(colour = "grey70", fill = "grey90", alpha = 0.6) +
    geom_vline(xintercept = quantile(data_df$z, probs=0.1), linetype="dashed") +
    geom_vline(xintercept = quantile(data_df$z, probs=0.9), linetype="dashed") +
    theme_void(base_size = 12)
  
  plot <- ggplot(plot_df, aes(x = .data[[var]], y = main_effect, group = TrajType, colour = TrajType, fill = TrajType)) +
    geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), alpha = 0.2) + 
    geom_line(alpha = 0.8, lwd = 0.6) +
    coord_cartesian(xlim = c(quantile(plot_df$var, probs=0.1),quantile(plot_df$var, probs=0.9)), ylim = c(0,1)) +
    geom_hline(yintercept= summary(analyses_df$increase)[[4]]/25, linetype='longdash', col = 'grey')+
    geom_hline(yintercept= summary(analyses_df$decrease)[[4]]/25, linetype='longdash', col = 'grey')+
    scale_fill_manual(values = c("#996633", "#669900")) + 
    scale_color_manual(values = c("#663300", "#61A36A")) +
    scale_x_continuous( labels = label_number()) +
    labs(x = x_label,
         y = "Probability of trajectory type") +
    theme_classic(base_size = 12) + theme(legend.position="none")
  
  x<- patchwork::wrap_elements(density) + 
    patchwork::wrap_elements(plot) +
    patchwork::plot_layout(ncol = 1, nrow = 2, widths = 4, heights = c(1, 4))
  x
  
}
burnedarea_plot <- plotting_function("trend_burnedarea", "trend in burned area (% per year)", burnedarea_plot_df )

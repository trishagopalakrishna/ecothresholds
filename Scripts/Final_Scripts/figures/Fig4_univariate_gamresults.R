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
increase_model_extractsmooths <- function (gam_model, outputsmooth_file_name_in_quotes){
  significant_onevariable_smooths <- gratia::smooth_estimates(gam_model,
                                                              select = c("s(anthropic_dist)", "s(trend_re)",
                                                                         "s(trend_annualtemp)",  "s(trend_heatwave)", 
                                                                         "s(trend_spei)", "s(trend_burnedarea)",
                                                                         "s(mean_re)", "s(mean_at)", "s(mean_spei)", 
                                                                         "s(mean_burnedarea)", "s(mean_savanna_percentage)",
                                                                         "s(mean_grass_percentage)", "s(mean_forest_percentage)"),
                                                              n=100, unconditional = TRUE) %>%
    mutate(intercept_no_uncertainty = gam_model$coefficients[1],
           intercept_add_estimate = .estimate + intercept_no_uncertainty,
           
           intercept_estimate_ci_low = intercept_add_estimate - .se*1.96,
           intercept_estimate_ci_high = intercept_add_estimate + .se*1.96,
           
           se_low_inv = inv.logit(intercept_estimate_ci_low),
           se_high_inv = inv.logit(intercept_estimate_ci_high),
           estimate_inv = inv.logit(intercept_add_estimate))
  significant_onevariable_smooths 
}
increase_significant_onevariablesmooths <- increase_model_extractsmooths(greening_gam)
increase_significant_onevariablesmooths <- increase_significant_onevariablesmooths %>% mutate(TrajType ="Greening")

decrease_model_extractsmooths <- function (gam_model, outputsmooth_file_name_in_quotes){
  significant_onevariable_smooths <- gratia::smooth_estimates(gam_model,
                                                              select = c("s(anthropic_dist)", "s(trend_re)",
                                                                         "s(trend_annualtemp)", "s(trend_spei)", "s(trend_burnedarea)",
                                                                         "s(mean_re)", "s(mean_at)", "s(mean_heatwave)",
                                                                         "s(mean_spei)", "s(mean_burnedarea)",
                                                                         "s(mean_grass_percentage)", "s(mean_forest_percentage)"),
                                                              n=100, unconditional = TRUE) %>%
    mutate(intercept_no_uncertainty = gam_model$coefficients[1],
           intercept_add_estimate = .estimate + intercept_no_uncertainty,
           
           intercept_estimate_ci_low = intercept_add_estimate - .se*1.96,
           intercept_estimate_ci_high = intercept_add_estimate + .se*1.96,
           
           se_low_inv = inv.logit(intercept_estimate_ci_low),
           se_high_inv = inv.logit(intercept_estimate_ci_high),
           estimate_inv = inv.logit(intercept_add_estimate))
  significant_onevariable_smooths
}
decrease_significant_onevariablesmooths <- decrease_model_extractsmooths(browning_gam)
decrease_significant_onevariablesmooths <- decrease_significant_onevariablesmooths %>% mutate(TrajType ="Browning")

#3. Plots
bothmodels_all_commonvar_labels <- read_csv(here("Outputs", "AnalysesResults", "GAMResults", "both_model_common_var_labels.csv"))
bothmodels_all_commonvar_list <- bothmodels_all_commonvar_labels$vars
bothmodels_all_commonvar_plots <- list(length = length(bothmodels_all_commonvar_list))

#3.(i) Plots for predictor variables that are significant in both browning and greening models
commonvar_plots <- purrr::map(bothmodels_all_commonvar_list, function(var){
  
  x_label <- bothmodels_all_commonvar_labels$label[bothmodels_all_commonvar_labels$vars == var]
  if (bothmodels_all_commonvar_labels$parse[bothmodels_all_commonvar_labels$vars == var]) {
    x_label <- parse(text = x_label)
  }
  
  data_df <- analyses_df_noNA_nooutliers %>%
    dplyr::select(all_of(var)) %>%
    dplyr::rename(z = 1)
  density <- ggplot(data_df, aes(x = z)) +
    geom_density(colour = "grey70", fill = "grey90", alpha = 0.6) +
    geom_vline(xintercept = quantile(data_df$z, probs=0.1), linetype="dashed") +
    geom_vline(xintercept = quantile(data_df$z, probs=0.9), linetype="dashed") +
    theme_void()
  
  increase_smooth_df <- increase_significant_onevariablesmooths  %>%
    dplyr::select(.smooth, .type, estimate_inv, se_low_inv, se_high_inv, all_of(var), TrajType) %>%
    rename(var0 = 6) %>%
    drop_na() %>% 
    filter(var0 > quantile(var0, probs=0.1) & var0 < quantile (var0, probs =0.9))
  
  decrease_smooth_df <- decrease_significant_onevariablesmooths  %>%
    dplyr::select(.smooth, .type, estimate_inv, se_low_inv, se_high_inv, all_of(var), TrajType) %>%
    rename(var0 = 6) %>%
    drop_na() %>%
    filter(var0 > quantile(var0, probs=0.1) & var0 < quantile (var0, probs =0.9))
  smooth_df <- bind_rows(increase_smooth_df, decrease_smooth_df)  
  
  plot <- ggplot(smooth_df, aes(x = var0, y = estimate_inv, group = TrajType, colour = TrajType, fill = TrajType)) +
    geom_ribbon(aes(ymin = se_low_inv, ymax = se_high_inv), alpha = 0.2) + 
    geom_line(alpha = 0.8, lwd = 0.6) + 
    coord_cartesian(xlim = c(quantile(smooth_df$var0, probs=0.1),quantile(smooth_df$var0, probs=0.9)), ylim = c(0,1)) +
    scale_fill_manual(values = c("#996633", "#669900"), name ="Trajectory Type") + 
    scale_color_manual(values = c("#663300", "#61A36A")) +
    scale_x_continuous( labels = label_number()) +
    labs(x = x_label,
         y = "Probability of trajectory type") +
    theme_classic() + theme(legend.position="none")
  
  x<-density + plot +
    patchwork::plot_layout(ncol = 1, nrow = 2, widths = 4, heights = c(1, 4))
})

#3.(ii) Plots of predictor variables that are significant only in one type of model
additional_plots <- function (var, model_smooth_df, x_label, fill_color, line_color){
  data_df <- analyses_df_noNA_nooutliers %>%
    dplyr::select(var) %>%
    dplyr::rename(z = 1)
  density <- ggplot(data_df, aes(x = z)) +
    geom_density(colour = "grey70", fill = "grey90", alpha = 0.6) +
    geom_vline(xintercept = quantile(data_df$z, probs=0.1), linetype="dashed") +
    geom_vline(xintercept = quantile(data_df$z, probs=0.9), linetype="dashed") +
    theme_void() 
  
  smooth_df <- model_smooth_df  %>%
    dplyr::select(.smooth, .type, estimate_inv, se_low_inv, se_high_inv, all_of(var), TrajType) %>%
    rename(var0 = 6) %>%
    drop_na() %>%
    filter(var0 > quantile(var0, probs=0.1) & var0 < quantile (var0, probs =0.9))
  
  plot <- ggplot(smooth_df, aes(x = var0, y = estimate_inv, group = TrajType, colour = TrajType, fill = TrajType)) +
    geom_ribbon(aes(ymin = se_low_inv, ymax = se_high_inv), alpha = 0.2) + 
    geom_line(alpha = 0.8, lwd = 0.6) + 
    coord_cartesian(xlim = c(quantile(smooth_df$var0, probs=0.1),quantile(smooth_df$var0, probs=0.9)), ylim = c(0,1)) +
    scale_fill_manual(values = fill_color, name ="Trajectory Type") + 
    scale_color_manual(values = line_color) +
    scale_x_continuous( labels = label_number()) +
    labs(x = x_label,
         y = "Probability of trajectory type") +
    theme_classic() + theme(legend.position="none")
  
  x <- density + plot +
    patchwork::plot_layout(ncol = 1, nrow = 2, widths = 4, heights = c(1, 4))
  x
}
mean_savanna_increase <- additional_plots("mean_savanna_percentage", increase_significant_onevariablesmooths, "mean savanna formation proportion (%)", "#669900", "#61A36A")
trend_heatwave_increase <- additional_plots("trend_heatwave", increase_significant_onevariablesmooths, "trend in heatwaves","#669900", "#61A36A")
mean_heatwave_decrease <- additional_plots("mean_heatwave", decrease_significant_onevariablesmooths, "mean number of heatwaves", "#996633", "#663300")

#3.(iii) Mean annual temperature predictor variable plot due to different x and y axis ranges
mean_annual_temp <- function (var, x_label){
  data_df <- analyses_df_noNA_nooutliers %>%
    dplyr::select(all_of(var)) %>%
    dplyr::rename(z = 1)
  density <- ggplot(data_df, aes(x = z)) +
    geom_density(colour = "grey70", fill = "grey90", alpha = 0.6) +
    geom_vline(xintercept = quantile(data_df$z, probs=0.1), linetype="dashed") +
    geom_vline(xintercept = quantile(data_df$z, probs=0.9), linetype="dashed") +
    theme_void()
  
  increase_smooth_df <- increase_significant_onevariablesmooths  %>%
    dplyr::select(.smooth, .type, estimate_inv, se_low_inv, se_high_inv, all_of(var), TrajType) %>%
    rename(var0 = 6) %>%
    drop_na() %>% 
    filter(var0 > quantile(var0, probs=0.1) & var0 < quantile (var0, probs =0.9))
  
  decrease_smooth_df <- decrease_significant_onevariablesmooths  %>%
    dplyr::select(.smooth, .type, estimate_inv, se_low_inv, se_high_inv, all_of(var), TrajType) %>%
    rename(var0 = 6) %>%
    drop_na() %>%
    filter(var0 > quantile(var0, probs=0.1) & var0 < quantile (var0, probs =0.9))
  smooth_df <- bind_rows(increase_smooth_df, decrease_smooth_df)  
  
  plot <- ggplot(smooth_df, aes(x = var0, y = estimate_inv, group = TrajType, colour = TrajType, fill = TrajType)) +
    geom_ribbon(aes(ymin = se_low_inv, ymax = se_high_inv), alpha = 0.2) + 
    geom_line(alpha = 0.8, lwd = 0.6) + 
    coord_cartesian(xlim = c(quantile(smooth_df$var0, probs=0.1),quantile(smooth_df$var0, probs=0.9)), ylim = c(0,1)) +
    scale_fill_manual(values = c("#996633", "#669900"), name ="Trajectory Type") + 
    scale_color_manual(values = c("#663300", "#61A36A")) +
    scale_x_continuous(breaks= c(295,296,297,298,299), labels = c(21.85, 22.85, 23.85, 24.85, 25.85))  +
    labs(x = x_label,
         y = "Probability of trajectory type") +
    theme_classic() + theme(legend.position="none")
  
  x<-density + plot +
    patchwork::plot_layout(ncol = 1, nrow = 2, widths = 4, heights = c(1, 4))
}
mean_at_plot <- mean_annual_temp ("mean_at", "mean annual temperature (degrees centigrade)")

#4. Combining all plots
plots <- ggpubr::ggarrange(plotlist = list(mean_savanna_increase, commonvar_plots[[10]], commonvar_plots[[9]], commonvar_plots[[1]],
                           mean_at_plot, commonvar_plots[[6]], mean_heatwave_decrease, commonvar_plots[[7]],
                           commonvar_plots[[5]], commonvar_plots[[2]], trend_heatwave_increase, commonvar_plots[[3]],
                           commonvar_plots[[8]], commonvar_plots[[4]]),
                           labels =c("a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)", "i)", "j)", "k)", "l)", "m)", "n)"),
                           ncol = 4, nrow =4,common.legend = F)
ggsave(here("Outputs", "TrendsResults", "Figures", "fig3.png"),
       plots,dpi = 700, height = 30, width = 30, units = "cm")
                                               
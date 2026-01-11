require(mgcv)
library(dplyr)
library(tidyverse)
library(tidyr)

browning_gam <- read_rds(here("Scripts", "Final_Scripts", "analyses" ,"onemodel_dec.rds")) #load browning gam model with your file path
analyses_df <-read_rds(here("Scripts", "Final_Scripts", "analyses" ,"analyses_df.rds")) #load data to generate new data grid


xgrid <- seq(min(analyses_df$trend_burnedarea),max(analyses_df$trend_burnedarea),len = 200)

zgrid <- c((mean(analyses_df$mean_savanna_percentage) - sd(analyses_df$mean_savanna_percentage)),
                 mean(analyses_df$mean_savanna_percentage), 
                 (mean(analyses_df$mean_savanna_percentage) + sd(analyses_df$mean_savanna_percentage))) 
newd <- data.frame( expand.grid(xgrid,zgrid) )
names(newd) <- c("trend_burnedarea", "mean_savanna_percentage") #match names in gam model

df_means<- analyses_df %>% dplyr::select(-c(x,y,"trend_burnedarea", "mean_savanna_percentage")) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm= T)))

x1 <- -43.2781
y1 <- -4.0213

newd <- newd%>% bind_cols(df_means) %>% mutate(x= x1, y=y1)


pred <- predict(browning_gam, newdata = newd, type="terms")
#removing lat, long information 
pred <- as.data.frame(pred)
xy_interaction_effects <- grep("x,y",names(pred),fixed = T)
pred_without_xy_interaction_effects <- (subset(pred, select = -xy_interaction_effects))

x1Eff_trial <- rowSums( pred_without_xy_interaction_effects[1:200, c(7,13, 22)] ) #1,13,18 are the columns associated with (anthropic_dist), (mean_savanna_percentage), (their interaction) from names(pred)
plot(xgrid,x1Eff_trial, type="l")
x2Eff_trial <- rowSums( pred_without_xy_interaction_effects[201:400, c(7,13,22)] )
lines(xgrid,x2Eff_trial, col = "red")
x3Eff_trial <- rowSums( pred_without_xy_interaction_effects[401:600, c(7,13,22)] )
lines(xgrid,x3Eff_trial, col = "blue")


# with uncertainty
X <- mgcv::predict.bam(browning_gam, newdata = newd, type ="lpmatrix")

###### THEO- THE GAM HERE SHOULD EXCLUDE LAT,LONG CORRECT? HAVE I DONE IT CORRECTLY IN LINE 55?
# beta uncertainty
n.sims <- 1000
b_sims <- rmvn(n.sims,coef(browning_gam),browning_gam$Vp) #if we continue keeping lat,long

#if we exclude lat,long from gam model 
xy_coef <- grep("x,y", names(coef(browning_gam)), fixed =T)

#selecting only required elements from variance-covariance matrix
#vp_matrix <- as.data.frame(browning_gam$Vp)
select_vp_matrix <- browning_gam$Vp[-xy_coef,-xy_coef]
b_sims_without_latlong <- rmvn(n.sims,coef(browning_gam)[-xy_coef],select_vp_matrix) 

# pick (anthropic_dist), (mean_savanna_percentage), (their interaction) effects only without lat,long
coef_index_trial1 <- grep("(trend_burnedarea)",names(coef(browning_gam)[-xy_coef]),fixed = T )
coef_index_trial2 <- grep("(mean_savanna_percentage)",names(coef(browning_gam)[-xy_coef]),fixed = T )
coef_index_trial3 <- grep("trend_burnedarea,mean_savanna_percentage",names(coef(browning_gam)[-xy_coef]),fixed = T )

coef_index<- c(coef_index_trial1, coef_index_trial2, coef_index_trial3)

# compute the spline
#x1Eff <- tcrossprod( X[,coef_index] , b_sims[,coef_index] )
X1Eff_without_latlong <- tcrossprod( X[,coef_index] , b_sims_without_latlong[,coef_index] )

######- PLOTS DO NOT WORK- THE "MEAN" PLOT IS THE SAME FOR ALL 3 VALUES OF MEAN_SAVANNA_PERCENTAGE WHICH IS NOT CORRECT AS SEEN ABOVE IN BLACK,RED,BLUE LINE PLOT
plot(xgrid, apply(X1Eff_without_latlong[1:200,],1,mean),type="l") #should be very similar to black line from line #34, but its not
lines(xgrid, apply(X1Eff_without_latlong[1:200,], 1, quantile,probs = 0.025), type = "l", lty = 3) 
lines(xgrid, apply(X1Eff_without_latlong[1:200,], 1, quantile,probs = 0.975), type = "l", lty = 3)

plot(xgrid, apply(X1Eff_without_latlong[201:400,],1,mean),type="l")#should be very similar to red line from line #36, but its not
lines(xgrid, apply(X1Eff_without_latlong[201:400,], 1, quantile,probs = 0.025), type = "l", lty = 3)
lines(xgrid, apply(X1Eff_without_latlong[201:400,], 1, quantile,probs = 0.975), type = "l", lty = 3)

plot(xgrid, apply(X1Eff_without_latlong[401:600,],1,mean),type="l")#should be very similar to blue line from line #38, but its not
lines(xgrid, apply(X1Eff_without_latlong[401:600,], 1, quantile,probs = 0.025), type = "l", lty = 3)
lines(xgrid, apply(X1Eff_without_latlong[401:600,], 1, quantile,probs = 0.975), type = "l", lty = 3)


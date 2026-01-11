require(mgcv)
library(dplyr)
library(tidyverse)
library(tidyr)

browning_gam <- read_rds(here("Scripts", "Final_Scripts", "analyses" ,"onemodel_dec.rds")) #load browning gam model with your file path
analyses_df <-read_rds(here("Scripts", "Final_Scripts", "analyses" ,"analyses_df.rds")) #load data to generate new data grid


xgrid <- seq(min(analyses_df$anthropic_dist),max(analyses_df$anthropic_dist),len = 200)

zgrid <- c((mean(analyses_df$mean_savanna_percentage) - sd(analyses_df$mean_savanna_percentage)),
                 mean(analyses_df$mean_savanna_percentage), 
                 (mean(analyses_df$mean_savanna_percentage) + sd(analyses_df$mean_savanna_percentage))) 
newd <- data.frame( expand.grid(xgrid,zgrid) )
names(newd) <- c("anthropic_dist", "mean_savanna_percentage") #match names in gam model

df_means<- analyses_df %>% dplyr::select(-c(x,y,"anthropic_dist", "mean_savanna_percentage")) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm= T)))

x1 <- -43.2781
y1 <- -4.0213

newd <- newd %>% bind_cols(df_means) %>% mutate(x= x1, y=y1)


pred <- predict(browning_gam, newdata = newd, type="terms")
#removing lat, long information 
pred <- as.data.frame(pred)
xy_interaction_effects <- grep("x,y",names(pred),fixed = T)
pred_without_xy_interaction_effects <- (subset(pred, select = -xy_interaction_effects))

########################### Theo- have I added the intercept and applied the inv.logit correctly?
x1Eff_trial <- rowSums( pred_without_xy_interaction_effects[1:200, c(1,13, 18)] ) #1,13,18 are the columns associated with (anthropic_dist), (mean_savanna_percentage), (their interaction) from names(pred)
x1Eff_intercept_invlogit <- inv.logit(x1Eff_trial+browning_gam$coefficients[1])
plot(xgrid,x1Eff_intercept_invlogit, type="l")

x2Eff_trial <- rowSums( pred_without_xy_interaction_effects[201:400, c(1,13, 18)] )
x2Eff_intercept_invlogit <- inv.logit(x2Eff_trial+browning_gam$coefficients[1])
lines(xgrid,x2Eff_intercept_invlogit, type = "l", col = "red")

x3Eff_trial <- rowSums( pred_without_xy_interaction_effects[401:600, c(1,13, 18)] )
x3Eff_intercept_invlogit <- inv.logit(x3Eff_trial +browning_gam$coefficients[1])
lines(xgrid, x3Eff_intercept_invlogit, col = "blue", type = "l")

# with uncertainty
X <- mgcv::predict.bam(browning_gam, newdata = newd, type ="lpmatrix")

n.sims <- 1000
b_sims <- rmvn(n.sims,coef(browning_gam),browning_gam$Vp) 

# pick (anthropic_dist), (mean_savanna_percentage), (their interaction) effects only without lat,long
coef_index_trial1 <- grep("s(anthropic_dist)",names(coef(browning_gam)),fixed = T )
coef_index_trial2 <- grep("s(mean_savanna_percentage)",names(coef(browning_gam)),fixed = T )
coef_index_trial3 <- grep("ti(anthropic_dist,mean_savanna_percentage)",names(coef(browning_gam)),fixed = T )

coef_index <- c(coef_index_trial1, coef_index_trial2, coef_index_trial3)

# compute the spline
X1Eff_without_latlong <- tcrossprod( X[,coef_index] , b_sims[,coef_index] )

mean_x1Eff <- apply(X1Eff_without_latlong[1:200,],1,mean)
invlogit_mean_x1Eff <- inv.logit(mean_x1Eff+browning_gam$coefficients[1])
plot(xgrid, invlogit_mean_x1Eff,type="l") 
low_ci_x1Eff <- apply(X1Eff_without_latlong[1:200,], 1, quantile,probs = 0.025)
invlogit_low_ci_x1Eff <- inv.logit(low_ci_x1Eff+browning_gam$coefficients[1])
lines(xgrid, invlogit_low_ci_x1Eff, type = "l", lty = 3) 
high_ci_x1Eff <- apply(X1Eff_without_latlong[1:200,], 1, quantile,probs = 0.975)
invlogit_high_ci_x1Eff <- inv.logit(high_ci_x1Eff+browning_gam$coefficients[1])
lines(xgrid, invlogit_high_ci_x1Eff, type = "l", lty = 3)

mean_x1Eff_2 <- apply(X1Eff_without_latlong[201:400,],1,mean)
invlogit_mean_x1Eff_2 <- inv.logit(mean_x1Eff_2+browning_gam$coefficients[1])
plot(xgrid, invlogit_mean_x1Eff_2,type="l")
low_ci_x1Eff_2 <- apply(X1Eff_without_latlong[201:400,], 1, quantile,probs = 0.025)
invlogit_low_ci_x1Eff_2 <- inv.logit(low_ci_x1Eff_2+browning_gam$coefficients[1])
lines(xgrid, invlogit_low_ci_x1Eff_2 , type = "l", lty = 3)
high_ci_x1Eff_2 <- apply(X1Eff_without_latlong[201:400,], 1, quantile,probs = 0.975)
invlogit_high_ci_x1Eff_2 <- inv.logit(high_ci_x1Eff_2+browning_gam$coefficients[1])
lines(xgrid, invlogit_high_ci_x1Eff_2, type = "l", lty = 3)

mean_x1Eff_3 <- apply(X1Eff_without_latlong[401:600,],1,mean)
invlogit_mean_x1Eff_3 <- inv.logit(mean_x1Eff_3+browning_gam$coefficients[1])
plot(xgrid, invlogit_mean_x1Eff_3, type="l")
low_ci_x1Eff_3 <- apply(X1Eff_without_latlong[401:600,], 1, quantile,probs = 0.025)
invlogit_low_ci_x1Eff_3 <- inv.logit(low_ci_x1Eff_3+browning_gam$coefficients[1])
lines(xgrid, invlogit_low_ci_x1Eff_3, type = "l", lty = 3)
high_ci_x1Eff_3 <- apply(X1Eff_without_latlong[401:600,], 1, quantile,probs = 0.975)
invlogit_high_ci_x1Eff_3 <- inv.logit(high_ci_x1Eff_3+browning_gam$coefficients[1])
lines(xgrid, invlogit_high_ci_x1Eff_3, type = "l", lty = 3)


library(tidyverse)
library(MuMIn)

#################################### Visual validation
class_trajectory_mod <- function (dataset = NULL, interval_size = 0.5) {
  dataset<- dataset %>% mutate(Months= 1:nrow(dataset))
  
  message("Trajectory fitting")
  null_mod<- lm(mean_valueint ~ 1, data = dataset) ##1. No intercept
  summary(null_mod) #significant pvalue next to intercept  means that the mean is significantly different from 0
  nmrs_null<- sqrt(sum(summary(null_mod)$residuals^2)/length(dataset$mean_valueint))/sd(dataset$mean_valueint)
  aic_null<- MuMIn::AICc(null_mod)
  
  linear_mod<- lm(mean_valueint ~ Months, data = dataset) ##2.linear
  summary(linear_mod)
  nmrs_lin<- sqrt(sum(summary(linear_mod)$residuals^2)/length(dataset$mean_valueint))/sd(dataset$mean_valueint)
  aic_lin<- MuMIn::AICc(linear_mod)
  
  orth_poly_mod<- lm(mean_valueint ~ poly(Months,2, raw=F), data = dataset) #raw=F means orthogonal polynomials are used
  summary(orth_poly_mod)
  nmrs_quad<- sqrt(sum(summary(orth_poly_mod)$residuals^2)/length(dataset$mean_valueint))/sd(dataset$mean_valueint)
  aic_quad<- MuMIn::AICc(orth_poly_mod)
  
  ## To get relevant values for quadratic output - Directly from Pellise et al., 2024 
  ## https://github.com/matpelissie/abrupt_shifts_ecological_timeseries_classification/blob/main/R/functions_trajclass.R
  
  # After getting Y = gamma*chi + delta*X' + epsilon with orthogonal polynomial,
  # we have to perform a variable change to obtain relevant values in the X interval 
  # for first_order_coefficient, second_order_coefficient and intercept, 
  # knowing that X'= alpha*X + beta and chi = eta*X'^2 + theta
  
  gammab  <-  orth_poly_mod$coefficients[3]
  delta  <-  orth_poly_mod$coefficients[2]
  epsilon  <-  orth_poly_mod$coefficients[1]
  
  alpha  <-  lm(orth_poly_mod$model[, 2][, 1] ~ dataset$Months)$coef[2]
  beta  <-  lm(orth_poly_mod$model[, 2][, 1] ~ dataset$Months)$coef[1]
  
  eta  <-  1/lm((orth_poly_mod$model[, 2][, 1])^2 ~
                  orth_poly_mod$model[, 2][, 2])$coef[2]
  theta  <-  (-lm((orth_poly_mod$model[, 2][, 1])^2 ~
                    orth_poly_mod$model[, 2][, 2])$coef[1])*eta
  
  Y2 <- dataset$mean_valueint*(max(dataset$Months)-min(dataset$Months))/(max(dataset$mean_valueint)-min(dataset$mean_valueint)) 
  
  # p2 and p3 are relevant when Y and X amplitudes are equivalent,in particular when 
  # studying scaled-to-1 indices, Y and X amplitudes may be very different, so we 
  # scaled the amplitudes to calculate p2 and p3
  
  polynomial_orthonormal_basis <- lm(Y2~poly(dataset$Months,2, raw=T))$coefficients
  
  message("Quadratic model output")
  classification <-
    data.frame(first_order_coefficient = (delta+2*beta*gammab*eta)*alpha,
               first_order_pvalue =
                 summary(orth_poly_mod)$coefficients[2, 4],
               second_order_coefficient = (alpha^2)*gammab*eta,
               second_order_pvalue =
                 summary(orth_poly_mod)$coefficients[3, 4],
               strd_error=summary(orth_poly_mod)$coefficients[2, 2],
               intercept = epsilon+beta*delta+(beta^2)*gammab*eta+gammab*theta,
               x_m = (dataset$Months[length(dataset$Months)]-dataset$Months[1])/2+dataset$Months[1],
               # points of interest:
               p1 = -(delta+2*beta*gammab*eta)/(2*alpha*gammab*eta),
               p2 = (-polynomial_orthonormal_basis[2]+1)/
                 (2*polynomial_orthonormal_basis[3]),
               p3 = (-polynomial_orthonormal_basis[2]-1)/
                 (2*polynomial_orthonormal_basis[3]),
               aic = aic_quad,
               nrmse = nmrs_quad,
               loc_brk = NA, 
               trend = NA,
               mag = NA,
               rel_chg = NA,
               SDbef = NA,
               SDaft = NA)
  
  message("Linear model output")
  classification[2,] <-
    data.frame(first_order_coefficient = delta*alpha,
               first_order_pvalue =
                 summary(orth_poly_mod)$coefficients[2, 4],
               second_order_coefficient = 0,
               second_order_pvalue =
                 summary(orth_poly_mod)$coefficients[3, 4],
               strd_error=summary(orth_poly_mod)$coefficients[2, 2],
               intercept = epsilon+delta*beta,
               x_m = (dataset$Months[length(dataset$Months)]-dataset$Months[1])/2+dataset$Months[1],
               p1 = NA,
               p2 = NA,
               p3 = NA,
               aic = aic_lin,
               nrmse = nmrs_lin,
               loc_brk = NA, 
               trend = NA,
               mag = NA,
               rel_chg = NA,
               SDbef = NA,
               SDaft = NA)
  
  message("No change model output")
  classification[3,] <-
    data.frame(first_order_coefficient = 0,
               first_order_pvalue =
                 summary(orth_poly_mod)$coefficients[2, 4],
               second_order_coefficient = 0,
               second_order_pvalue =
                 summary(orth_poly_mod)$coefficients[3, 4],
               strd_error=summary(orth_poly_mod)$coefficients[2, 2],
               intercept = null_mod$coefficients,
               x_m = (dataset$Months[length(dataset$Months)]-dataset$Months[1])/2+dataset$Months[1],
               p1 = NA,
               p2 = NA,
               p3 = NA,
               aic = aic_null,
               nrmse = nmrs_null,
               loc_brk = NA, 
               trend = NA,
               mag = NA,
               rel_chg = NA,
               SDbef = NA,
               SDaft = NA)
  
  
  message("Classification of each model fitted")
  for (i in 1:3){
    
    # Compute the derivative at xm-delta and at xm + delta with delta being
    # half of the input interval size (its 25% of the time interval which is set as 0.5)
    derivative <-
      2*(classification$x_m[i] - (dataset$Months[length(dataset$Months)]-dataset$Months[1])*(interval_size/2))*  
      classification$second_order_coefficient[i] +
      classification$first_order_coefficient[i]
    derivative2 <-
      2*(classification$x_m[i] + (dataset$Months[length(dataset$Months)]-dataset$Months[1])*(interval_size/2))*
      classification$second_order_coefficient[i] +
      classification$first_order_coefficient[i]
    
    
    if(sign(derivative) != sign(derivative2) | i==3){
      # non consistent direction around x_m
      classification$derivative[i]  <-  NA
      classification$intercept_derivative[i]  <-  NA
      
    } else {
      # consistent direction around x_m
      classification$derivative[i] <- mean(c(derivative, derivative2))
      classification$intercept_derivative[i] <-
        (classification$second_order_coefficient[i]*classification$x_m[i]^2+
           classification$first_order_coefficient[i]*classification$x_m[i]+
           classification$intercept[i]) -
        classification$x_m[i]*classification$derivative[i]
    }
    
    # Compute the derivative of the curvature function to get acceleration:
    classification$derivated_curvature[i] <-
      -12*(classification$second_order_coefficient[i]^2)*
      (2*classification$second_order_coefficient[i]*classification$x_m[i]+
         classification$first_order_coefficient[i])*
      (classification$second_order_coefficient[i]/
         abs(classification$second_order_coefficient[i]))/
      ((1+(2*classification$second_order_coefficient[i]*classification$x_m[i]+
             classification$first_order_coefficient[i])^2)^(2.5))
    
    # Keep derivated curvature even if not significant for polynomial fit:
    if(classification$second_order_pvalue[i]>0.05 & i != 1){
      classification$derivated_curvature[i] <- NA
    }
    
    # Classify the direction:
    classification$direction[i] <- NA
    classification$direction[i][which(
      classification$derivative[i] > 0)] <- "increase"
    classification$direction[i][which(
      classification$derivative[i] < 0)] <- "decrease"
    classification$direction[i][which(
      is.na(classification$derivative[i]))] <- "stable"
    classification$direction[i][which(
      as.numeric(classification$first_order_pvalue[i])>0.05 &
        as.numeric(classification$second_order_pvalue[i])>0.05)] <- "stable"
    
    # Classify the acceleration:
    classification$acceleration[i] <- NA
    classification$acceleration[i][which(
      classification$derivated_curvature[i] < 0)] <- "accelerated"
    classification$acceleration[i][which(
      classification$derivated_curvature[i] > 0)] <- "decelerated"
    classification$acceleration[i][which(
      classification$direction[i] == "stable" &
        classification$second_order_coefficient[i] < 0)] <- "concave"
    classification$acceleration[i][which(
      classification$direction[i] == "stable" &
        classification$second_order_coefficient[i] > 0)] <- "convex"
    classification$acceleration[i][which(
      is.na(classification$derivated_curvature[i]))] <- "constant"
    
    # Give the final classification combining direction and acceleration:
    classification$shape_class[i] <- paste(classification$direction[i],
                                           classification$acceleration[i],
                                           sep="_")
  }
  
  message ("Abrupt breakpoint analyses")
  chng_fit<- chngpt::chngptm(formula.1=mean_valueint~1,
                             formula.2=~Months,
                             family="gaussian", data=dataset,
                             type="step",
                             var.type="bootstrap", weights=NULL)
  
  pred_chg <- data.frame(timestep = dataset$Months,
                         bp = chng_fit$best.fit$fitted.values)
  
  chng_nrmse <- sqrt(sum(residuals(chng_fit)^2)/length(dataset$mean_valueint))/sd(dataset$mean_valueint)
  
  classification[4, 13] <- chng_fit$chngpt
  classification[4, 11] <- MuMIn::AICc(chng_fit)
  classification[4, 12] <- chng_nrmse
  classification[4, 14] <- ifelse(pred_chg$bp[1] >
                                    pred_chg$bp[length(pred_chg$bp)],
                                  "decrease", "increase")
  classification[4, 15] <-  pred_chg$bp[length(pred_chg$bp)] - pred_chg$bp[1]
  classification[4, 16] <-  (pred_chg$bp[length(pred_chg$bp)] - pred_chg$bp[1]) /
    max(abs(pred_chg$bp[length(pred_chg$bp)]),abs(pred_chg$bp[1]))
  classification[4, 17] <-  dataset %>%
    dplyr::filter(Months<=chng_fit$chngpt) %>%
    dplyr::pull(mean_valueint) %>%
    sd()
  classification[4, 18] <- dataset %>%
    dplyr::filter(Months>=chng_fit$chngpt) %>%
    dplyr::pull(mean_valueint) %>%
    sd()
  
  row.names(classification) <- c("Y_pol","Y_lin", "Y_nch", "Y_abpt")
  
  # classification$cell <- unique(dataset$cell)
  classification$x <- unique(dataset$x)
  classification$y <- unique(dataset$y)
  
  return(classification)
}
df_prep<- function (df){
  df<- df %>% separate_wider_delim(cols = name, delim = "_", names = c("Year","Month"))
  df<- df %>% group_by(Year) %>% summarise(trend = mean(trend, na.rm= TRUE))
  df
}
model_levels <- c("Null", "Lin", "Step", "Quad")

model_select<- function(models_pixel_df){
  message("Adding model name")
  models_pixel_df <- models_pixel_df %>% mutate(model_order=ordered(c("Quad", "Lin", "Null", "Step"),
                                                                    levels = model_levels))
  
  message("Least AIC model selection")
  models_pixel_df <- models_pixel_df %>%
    mutate(aic_diff= abs(aic) - min(abs(aic))) #wrt model with least AIC
  
  message("Conditional model selection")
  condition_less2<- models_pixel_df %>% filter(aic_diff<=2)
  
}

finaldf<- read.csv(here("Outputs", "TrajectoryPlotting", "sampled_points_timeseries_df.csv"))

sample<- read.csv(here("Outputs", "STL", "trialdata.csv"))
x_point<- sample %>% filter( cell== 4343430)
x_point<- df_prep(x_point)
traj_x_point <- class_trajectory_mod(x_point)
aic_condition<- model_select(traj_x_point)
remove(x_point, traj_x_point, aic_condition)

x_point<- sample %>% filter( cell== 4434277)
x_point<- df_prep(x_point)
traj_x_point <- class_trajectory_mod(x_point)
aic_condition<- model_select(traj_x_point)
remove(x_point, traj_x_point, aic_condition)

x_point<- sample %>% filter( cell== 4508588)
x_point<- df_prep(x_point)
traj_x_point <- class_trajectory_mod(x_point)
aic_condition<- model_select(traj_x_point)
remove(x_point, traj_x_point, aic_condition)

sample2<- read.csv(here("Outputs", "STL", "trialdata2.csv"))
x_point<- sample2 %>% filter( cell== 2542759)
x_point<- df_prep(x_point)
traj_x_point <- class_trajectory_mod(x_point)
aic_condition<- model_select(traj_x_point)
remove(x_point, traj_x_point, aic_condition)

sample3<- read.csv(here("Outputs", "STL", "trialdata3.csv"))
x_point<- sample3 %>% filter( cell== 853842)
x_point<- df_prep(x_point)
traj_x_point <- class_trajectory_mod(x_point)
aic_condition<- model_select(traj_x_point)
remove(x_point, traj_x_point, aic_condition)

sample4<- read.csv(here("Outputs", "STL", "trialdata4.csv"))
x_point<- sample4 %>% filter( cell== 1761954)
x_point<- df_prep(x_point)
traj_x_point <- class_trajectory_mod(x_point)
aic_condition<- model_select(traj_x_point)
remove(x_point, traj_x_point, aic_condition)

x_point<- sample4 %>% filter( cell== 2028251)
x_point<- df_prep(x_point)
traj_x_point <- class_trajectory_mod(x_point)
aic_condition<- model_select(traj_x_point)
remove(x_point, traj_x_point, aic_condition)


#################################### Checking classification funciton against Pelisse et al
################################### https://github.com/matpelissie/abrupt_shifts_ecological_timeseries_classification/blob/main/R/functions_trajclass.R

class_trajectory_mod <- function (dataset = NULL, interval_size = 0.5) {
  dataset<- dataset %>% rename("Y" = "trend")
  dataset<- dataset %>% mutate(X= 1:nrow(dataset))
  
  message("Trajectory fitting")
  null_mod<- lm(Y ~ 1, data = dataset) ##1. No intercept
  summary(null_mod) #significant pvalue neYt to intercept  means that the mean is significantly different from 0
  nmrs_null<- sqrt(sum(summary(null_mod)$residuals^2)/length(dataset$Y))/sd(dataset$Y)
  aic_null<- MuMIn::AICc(null_mod)
  
  linear_mod<- lm(Y ~ X, data = dataset) ##2.linear
  summary(linear_mod)
  nmrs_lin<- sqrt(sum(summary(linear_mod)$residuals^2)/length(dataset$Y))/sd(dataset$Y)
  aic_lin<- MuMIn::AICc(linear_mod)
  
  orth_poly_mod<- lm(Y ~ poly(X,2, raw=F), data = dataset) #raw=F means orthogonal polynomials are used
  summary(orth_poly_mod)
  nmrs_quad<- sqrt(sum(summary(orth_poly_mod)$residuals^2)/length(dataset$Y))/sd(dataset$Y)
  aic_quad<- MuMIn::AICc(orth_poly_mod)
  
  ## To get relevant values for quadratic output - Directly from Pellise et al., 2024 
  ## https://github.com/matpelissie/abrupt_shifts_ecological_timeseries_classification/blob/main/R/functions_trajclass.R
  
  # After getting Y = gamma*chi + delta*X' + epsilon with orthogonal polynomial,
  # we have to perform a variable change to obtain relevant values in the X interval 
  # for first_order_coefficient, second_order_coefficient and intercept, 
  # knowing that X'= alpha*X + beta and chi = eta*X'^2 + theta
  
  gammab  <-  orth_poly_mod$coefficients[3]
  delta  <-  orth_poly_mod$coefficients[2]
  epsilon  <-  orth_poly_mod$coefficients[1]
  
  alpha  <-  lm(orth_poly_mod$model[, 2][, 1] ~ dataset$X)$coef[2]
  beta  <-  lm(orth_poly_mod$model[, 2][, 1] ~ dataset$X)$coef[1]
  
  eta  <-  1/lm((orth_poly_mod$model[, 2][, 1])^2 ~
                  orth_poly_mod$model[, 2][, 2])$coef[2]
  theta  <-  (-lm((orth_poly_mod$model[, 2][, 1])^2 ~
                    orth_poly_mod$model[, 2][, 2])$coef[1])*eta
  
  Y2 <- dataset$Y*(max(dataset$X)-min(dataset$X))/(max(dataset$Y)-min(dataset$Y)) 
  
  # p2 and p3 are relevant when Y and X amplitudes are equivalent,in particular when 
  # studying scaled-to-1 indices, Y and X amplitudes may be very different, so we 
  # scaled the amplitudes to calculate p2 and p3
  
  polynomial_orthonormal_basis <- lm(Y2~poly(dataset$X,2, raw=T))$coefficients
  
  message("Quadratic model output")
  classification <-
    data.frame(first_order_coefficient = (delta+2*beta*gammab*eta)*alpha,
               first_order_pvalue =
                 summary(orth_poly_mod)$coefficients[2,4],
               second_order_coefficient = (alpha^2)*gammab*eta,
               second_order_pvalue =
                 summary(orth_poly_mod)$coefficients[3,4],
               strd_error=summary(orth_poly_mod)$coefficients[2,2],
               intercept = epsilon+beta*delta+(beta^2)*gammab*eta+gammab*theta,
               x_m = (dataset$X[length(dataset$X)]-dataset$X[1])/2+dataset$X[1],
               # points of interest:
               p1 = -(delta+2*beta*gammab*eta)/(2*alpha*gammab*eta),
               p2 = (-polynomial_orthonormal_basis[2]+1)/
                 (2*polynomial_orthonormal_basis[3]),
               p3 = (-polynomial_orthonormal_basis[2]-1)/
                 (2*polynomial_orthonormal_basis[3]),
               aic = aic_quad,
               nrmse = nmrs_quad,
               loc_brk = NA, 
               trend = NA,
               mag = NA,
               rel_chg = NA,
               SDbef = NA,
               SDaft = NA)
  
  message("Linear model output")
  classification[2,] <-
    data.frame(first_order_coefficient = delta*alpha,
               first_order_pvalue =
                 summary(orth_poly_mod)$coefficients[2, 4],
               second_order_coefficient = 0,
               second_order_pvalue =
                 summary(orth_poly_mod)$coefficients[3, 4],
               strd_error=summary(orth_poly_mod)$coefficients[2, 2],
               intercept = epsilon+delta*beta,
               x_m = (dataset$X[length(dataset$X)]-dataset$X[1])/2+dataset$X[1],
               p1 = NA,
               p2 = NA,
               p3 = NA,
               aic = aic_lin,
               nrmse = nmrs_lin,
               loc_brk = NA, 
               trend = NA,
               mag = NA,
               rel_chg = NA,
               SDbef = NA,
               SDaft = NA)
  
  message("No change model output")
  classification[3,] <-
    data.frame(first_order_coefficient = 0,
               first_order_pvalue =
                 summary(orth_poly_mod)$coefficients[2, 4],
               second_order_coefficient = 0,
               second_order_pvalue =
                 summary(orth_poly_mod)$coefficients[3, 4],
               strd_error=summary(orth_poly_mod)$coefficients[2, 2],
               intercept = null_mod$coefficients,
               x_m = (dataset$X[length(dataset$X)]-dataset$X[1])/2+dataset$X[1],
               p1 = NA,
               p2 = NA,
               p3 = NA,
               aic = aic_null,
               nrmse = nmrs_null,
               loc_brk = NA, 
               trend = NA,
               mag = NA,
               rel_chg = NA,
               SDbef = NA,
               SDaft = NA)
  
  message("Classification of each model fitted")
  for (i in 1:3){
    
    # Compute the derivative at xm-delta and at xm + delta with delta being
    # half of the input interval size (its 25% of the time interval which is set as 0.5)
    derivative <-
      2*(classification$x_m[i] - (dataset$X[length(dataset$X)]-dataset$X[1])*(interval_size/2))*  
      classification$second_order_coefficient[i] +
      classification$first_order_coefficient[i]
    derivative2 <-
      2*(classification$x_m[i] + (dataset$X[length(dataset$X)]-dataset$X[1])*(interval_size/2))*
      classification$second_order_coefficient[i] +
      classification$first_order_coefficient[i]
    
    
    if(sign(derivative) != sign(derivative2) | i==3){
      # non consistent direction around x_m
      classification$derivative[i]  <-  NA
      classification$intercept_derivative[i]  <-  NA
      
    } else {
      # consistent direction around x_m
      classification$derivative[i] <- mean(c(derivative, derivative2))
      classification$intercept_derivative[i] <-
        (classification$second_order_coefficient[i]*classification$x_m[i]^2+
           classification$first_order_coefficient[i]*classification$x_m[i]+
           classification$intercept[i]) -
        classification$x_m[i]*classification$derivative[i]
    }
    
    # Compute the derivative of the curvature function to get acceleration:
    classification$derivated_curvature[i] <-
      -12*(classification$second_order_coefficient[i]^2)*
      (2*classification$second_order_coefficient[i]*classification$x_m[i]+
         classification$first_order_coefficient[i])*
      (classification$second_order_coefficient[i]/
         abs(classification$second_order_coefficient[i]))/
      ((1+(2*classification$second_order_coefficient[i]*classification$x_m[i]+
             classification$first_order_coefficient[i])^2)^(2.5))
    
    # Keep derivated curvature even if not significant for polynomial fit:
    if(classification$second_order_pvalue[i]>0.05 & i != 1){
      classification$derivated_curvature[i] <- NA
    }
    
    # Classify the direction:
    classification$direction[i] <- NA
    classification$direction[i][which(
      classification$derivative[i] > 0)] <- "increase"
    classification$direction[i][which(
      classification$derivative[i] < 0)] <- "decrease"
    classification$direction[i][which(
      is.na(classification$derivative[i]))] <- "stable"
    classification$direction[i][which(
      as.numeric(classification$first_order_pvalue[i])>0.05 &
        as.numeric(classification$second_order_pvalue[i])>0.05)] <- "stable"
    
    # Classify the acceleration:
    classification$acceleration[i] <- NA
    classification$acceleration[i][which(
      classification$derivated_curvature[i] < 0)] <- "accelerated"
    classification$acceleration[i][which(
      classification$derivated_curvature[i] > 0)] <- "decelerated"
    classification$acceleration[i][which(
      classification$direction[i] == "stable" &
        classification$second_order_coefficient[i] < 0)] <- "concave"
    classification$acceleration[i][which(
      classification$direction[i] == "stable" &
        classification$second_order_coefficient[i] > 0)] <- "convex"
    classification$acceleration[i][which(
      is.na(classification$derivated_curvature[i]))] <- "constant"
    
    # Give the final classification combining direction and acceleration:
    classification$shape_class[i] <- paste(classification$direction[i],
                                           classification$acceleration[i],
                                           sep="_")
  }
  
  message ("Abrupt breakpoint analyses")
  chng_fit<- chngpt::chngptm(formula.1=Y~1,
                             formula.2=~X,
                             family="gaussian", data=dataset,
                             type="step",
                             var.type="bootstrap", weights=NULL)
  
  pred_chg <- data.frame(timestep = dataset$X,
                         bp = chng_fit$best.fit$fitted.values)
  
  chng_nrmse <- sqrt(sum(residuals(chng_fit)^2)/length(dataset$Y))/sd(dataset$Y)
  
  classification[4, 13] <- chng_fit$chngpt
  classification[4, 11] <- MuMIn::AICc(chng_fit)
  classification[4, 12] <- chng_nrmse
  classification[4, 14] <- ifelse(pred_chg$bp[1] >
                                    pred_chg$bp[length(pred_chg$bp)],
                                  "decrease", "increase")
  classification[4, 15] <-  pred_chg$bp[length(pred_chg$bp)] - pred_chg$bp[1]
  classification[4, 16] <-  (pred_chg$bp[length(pred_chg$bp)] - pred_chg$bp[1]) /
    max(abs(pred_chg$bp[length(pred_chg$bp)]),abs(pred_chg$bp[1]))
  classification[4, 17] <-  dataset %>%
    dplyr::filter(Y<=chng_fit$chngpt) %>%
    dplyr::pull(X) %>%
    sd()
  classification[4, 18] <- dataset %>%
    dplyr::filter(X>=chng_fit$chngpt) %>%
    dplyr::pull(X) %>%
    sd()
  
  row.names(classification) <- c("Y_pol","Y_lin", "Y_nch", "Y_abpt")
  
  # classification$cell <- unique(dataset$cell)
  classification$x <- unique(dataset$x)
  classification$y <- unique(dataset$y)
  
  return(classification)
} #Modified
model_select<- function(models_pixel_df){ #modified
  message("Adding model name")
  models_pixel_df <- models_pixel_df %>% mutate(model_order=ordered(c("Quad", "Lin", "Null", "Step"),
                                                                    levels = model_levels))
  
  message("Least AIC model selection")
  models_pixel_df <- models_pixel_df %>%
    mutate(aic_diff= aic - min(aic)) #wrt model with least AIC
  
  message("Conditional model selection")
  condition_less2<- models_pixel_df %>% filter(aic_diff<=2)
  
  if (dim(condition_less2)[1]==1){
    models_pixel_df<- condition_less2 
  } else{
    models_pixel_df<- condition_less2 %>% 
      filter(model_order==min(model_order))
  }
  models_pixel_df 
}


sample<- read.csv(here("Outputs", "STL", "trialdata.csv"))
x_point<- sample %>% filter( cell== 4343430)
x_point<- df_prep(x_point)
dataset<- x_point
x_traj<- class_trajectory_mod(dataset)
models_pixel_df<- x_traj
remove(x_point, dataset, x_traj, models_pixel_df, condition_less2)

x_point<- sample %>% filter( cell== 4434277)
x_point<- df_prep(x_point)
dataset<- x_point
x_traj<- class_trajectory_mod(dataset)
models_pixel_df<- x_traj
remove(x_point, dataset, x_traj, models_pixel_df, condition_less2)

x_point<- sample %>% filter( cell== 4508588)
x_point<- df_prep(x_point)
dataset<- x_point
x_traj<- class_trajectory_mod(dataset)
models_pixel_df <- x_traj
remove(x_point, dataset, x_traj, models_pixel_df, condition_less2)
remove(sample)

sample2<- read.csv(here("Outputs", "STL", "trialdata2.csv"))
x_point<- sample2 %>% filter( cell== 2542759)
x_point<- df_prep(x_point)
dataset<- x_point
x_traj<- class_trajectory_mod(dataset)
models_pixel_df <- x_traj
remove(x_point, dataset, x_traj, models_pixel_df, condition_less2)
remove(sample2)

sample3<- read.csv(here("Outputs", "STL", "trialdata3.csv"))
x_point<- sample3 %>% filter( cell== 853842)
x_point<- df_prep(x_point)
dataset<- x_point
x_traj<- class_trajectory_mod(dataset)
models_pixel_df <- x_traj
remove(x_point, dataset, x_traj, models_pixel_df, condition_less2)
remove(sample3)


sample4<- read.csv(here("Outputs", "STL", "trialdata4.csv"))
x_point<- sample4 %>% filter( cell== 1761954)
x_point<- df_prep(x_point)
dataset<- x_point
x_traj<- class_trajectory_mod(dataset)
models_pixel_df <- x_traj
remove(x_point, dataset, x_traj, models_pixel_df, condition_less2)


x_point<- sample4 %>% filter( cell== 2028251)
x_point<- df_prep(x_point)
dataset<- x_point
x_traj<- class_trajectory_mod(dataset)
models_pixel_df <- x_traj
remove(x_point, dataset, x_traj, models_pixel_df, condition_less2)
remove(sample4) #After modified trajectory and model selection function
# quadratic trajectory is selected, which I guess I can see. Next lowest
#AICdiff is Step decrease, which is what I thought is the trajectory
  

library(tidyverse)
library(foreach)
library(doParallel)
library(MuMIn)

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


#Data input
x<-read_rds("/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/anisoEVI/STL_Decomposition/swindow_11/stl_southern1_10.rds") #changed this manually

chunk2 <- 240*100
n2 <- nrow(x)
r2  <- rep(1:ceiling(n2/chunk2),each=chunk2)[1:(n2)]
split_pivot_df <- split(x, r2)
remove(chunk2, n2,r2)
gc()
message("Big rds file split")

numCores<- 16
library(foreach)
library(doParallel)

my.cluster<- parallel::makeCluster(
  numCores,
  type = "FORK",
  outfile = ""
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered() #should be TRUE
foreach::getDoParWorkers() #should give numCores


##Converting above to nested parallelized foreach loop
Sys.time();y<- foreach(
  i= 1:length(split_pivot_df),
  .combine= "rbind") %dopar%
  (split_pivot_df[[i]] %>%
     separate_wider_delim(name, "_", names=c("Year", "Month")) %>%
     group_by(cell,x,y,Year) %>% 
     summarise(mean_valueint=mean(value_int, na.rm=T), 
               mean_trend= mean(trend, na.rm=T)) %>%
     group_by(cell) %>%
     nest() %>%
     mutate(classification_data = purrr::map(data, class_trajectory_mod)) %>%
     select(-data) %>%
     unnest(classification_data)
  ); Sys.time()

write_rds(y, "/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/anisoEVI/TrajectoryShapes/annual/no_detrending/annualmean_trajshape_southern1_10.rds")#manually changed

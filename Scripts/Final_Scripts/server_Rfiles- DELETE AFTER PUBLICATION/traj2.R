library(tidyverse)
library(foreach)
library(doParallel)
library(MuMIn)
library(chngpt)


trajectory_classification <- function (dataset = NULL, value_int_or_trend, interval_size = 0.5) {
  dataset<- dataset %>% mutate(X= 1:nrow(dataset))
  #dataset<- dataset %>% rename("Y" = value_int_or_trend)
  dataset <- dataset %>% mutate(Y= dataset[[value_int_or_trend]])
  
  message("Trajectory fitting")
  null_mod<- lm(Y ~ 1, data = dataset) ##1. No intercept
  summary(null_mod) #significant pvalue next to intercept  means that the mean is significantly different from 0
  nmrs_null<- sqrt(sum(summary(null_mod)$residuals^2)/length(dataset$Y))/sd(dataset$Y)
  aic_null<- MuMIn::AICc(null_mod)
  
  linear_mod<- lm(Y ~ X, data = dataset) ##2.linear
  summary(linear_mod)
  nmrs_lin<- sqrt(sum(summary(linear_mod)$residuals^2)/length(dataset$Y))/sd(dataset$Y)
  aic_lin<- MuMIn::AICc(linear_mod)
  
  orth_poly_mod<- lm(Y ~ poly(X, 2, raw=F), data = dataset) #raw=F means orthogonal polynomials are used
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
                 summary(orth_poly_mod)$coefficients[2, 4],
               second_order_coefficient = (alpha^2)*gammab*eta,
               second_order_pvalue =
                 summary(orth_poly_mod)$coefficients[3, 4],
               strd_error=summary(orth_poly_mod)$coefficients[2, 2],
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
  chng_fit<- chngpt::chngptm(formula.1=Y ~ 1,
                             formula.2= ~X,
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
    dplyr::filter(X<=chng_fit$chngpt) %>%
    dplyr::pull(Y) %>%
    sd()
  classification[4, 18] <- dataset %>%
    dplyr::filter(X>=chng_fit$chngpt) %>%
    dplyr::pull(Y) %>%
    sd()
  
  row.names(classification) <- c("Y_pol","Y_lin", "Y_nch", "Y_abpt")
  
  # classification$cell <- unique(dataset$cell)
  classification$x <- unique(dataset$x)
  classification$y <- unique(dataset$y)
  
  return(classification)
}

numCores<- 20

my.cluster<- parallel::makeCluster(
  numCores,
  type = "FORK",
  outfile = ""
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered() #should be TRUE
foreach::getDoParWorkers() #should give numCores

index_file_path <- "/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/"

trajectory_function <- function (index_name, swinfolder_input, climate_zone, value_int_or_trend, swinfolder_output){
  file_list<- list.files (path = paste0(index_file_path, index_name,"/", "STL_Decomposition/", swinfolder_input, "/"), pattern =paste0("*",climate_zone,"*"), all.files = T, full.names = T)
  file_list <- gtools::mixedsort(file_list)
  message ("reading file paths complete")
  
  y<- foreach(
    i = 1:length(file_list)) %dopar% {
      read_rds(file_list[[i]]) %>% 
        separate_wider_delim(name, "_", names=c("Year", "Month")) %>% 
        group_by(cell) %>%
        nest() %>%
        mutate(classification_data = purrr::map(data, trajectory_classification, value_int_or_trend = value_int_or_trend)) %>%
        select(-data) %>%
        unnest(classification_data)
    }
  compiled_df <- bind_rows(y)
  write_rds(compiled_df, paste0(index_file_path, index_name, "/", "TrajectoryShapes/monthly/",swinfolder_output, "/", "trajectory_", climate_zone, ".rds"))
  message ("file written to disk")
}

#NDVI
#Sys.time(); trajectory_function("NDVI", "swindow_7", "central", "trend", "swindow_7"); Sys.time()
#Sys.time(); trajectory_function("NDVI", "swindow_7", "eastern", "trend", "swindow_7"); Sys.time()
#Sys.time(); trajectory_function("NDVI", "swindow_7", "southern", "trend", "swindow_7"); Sys.time()

#Sys.time(); trajectory_function("NDVI", "swindow_11", "central", "trend", "swindow_11"); Sys.time()
#Sys.time(); trajectory_function("NDVI", "swindow_11", "eastern", "trend", "swindow_11"); Sys.time()
#Sys.time(); trajectory_function("NDVI", "swindow_11", "southern", "trend", "swindow_11"); Sys.time()

#Sys.time(); trajectory_function("NDVI", "swindow_periodic", "central", "trend", "swindow_periodic"); Sys.time()
#Sys.time(); trajectory_function("NDVI", "swindow_periodic", "eastern", "trend","swindow_periodic"); Sys.time()
#Sys.time(); trajectory_function("NDVI", "swindow_periodic", "southern", "trend", "swindow_periodic"); Sys.time()

#Sys.time(); trajectory_function("NDVI", "swindow_periodic", "central", "value_int", "value_int"); Sys.time()
#Sys.time(); trajectory_function("NDVI", "swindow_periodic", "eastern", "value_int", "value_int"); Sys.time()
#Sys.time(); trajectory_function("NDVI", "swindow_periodic", "southern", "value_int", "value_int"); Sys.time()

#EVI
#Sys.time(); trajectory_function("EVI", "swindow_7", "central", "trend", "swindow_7"); Sys.time()
#Sys.time(); trajectory_function("EVI", "swindow_7", "eastern", "trend", "swindow_7"); Sys.time()
#Sys.time(); trajectory_function("EVI", "swindow_7", "southern", "trend","swindow_7"); Sys.time()

#Sys.time(); trajectory_function("EVI", "swindow_11", "central", "trend", "swindow_11"); Sys.time()
#Sys.time(); trajectory_function("EVI", "swindow_11", "eastern", "trend","swindow_11"); Sys.time()
#Sys.time(); trajectory_function("EVI", "swindow_11", "southern", "trend", "swindow_11"); Sys.time()

#Sys.time(); trajectory_function("EVI", "swindow_periodic", "central", "trend", "swindow_periodic"); Sys.time()
#Sys.time(); trajectory_function("EVI", "swindow_periodic", "eastern", "trend", "swindow_periodic"); Sys.time()
#Sys.time(); trajectory_function("EVI", "swindow_periodic", "southern", "trend", "swindow_periodic"); Sys.time()

#Sys.time(); trajectory_function("EVI", "swindow_periodic", "central", "value_int", "value_int"); Sys.time()
#Sys.time(); trajectory_function("EVI", "swindow_periodic", "eastern", "value_int", "value_int"); Sys.time()
#Sys.time(); trajectory_function("EVI", "swindow_periodic", "southern", "value_int", "value_int"); Sys.time()

#anisoEVI
Sys.time(); trajectory_function("anisoEVI", "swindow_7", "central", "trend", "swindow_7"); Sys.time()
Sys.time(); trajectory_function("anisoEVI", "swindow_7", "eastern", "trend", "swindow_7"); Sys.time()
Sys.time(); trajectory_function("anisoEVI", "swindow_7", "southern", "trend", "swindow_7"); Sys.time()

Sys.time(); trajectory_function("anisoEVI", "swindow_11", "central", "trend","swindow_11"); Sys.time()
Sys.time(); trajectory_function("anisoEVI", "swindow_11", "eastern", "trend", "swindow_11"); Sys.time()
Sys.time(); trajectory_function("anisoEVI", "swindow_11", "southern", "trend", "swindow_11"); Sys.time()

Sys.time(); trajectory_function("anisoEVI", "swindow_periodic", "central", "trend", "swindow_periodic"); Sys.time()
Sys.time(); trajectory_function("anisoEVI", "swindow_periodic", "eastern", "trend", "swindow_periodic"); Sys.time()
Sys.time(); trajectory_function("anisoEVI", "swindow_periodic", "southern", "trend", "swindow_periodic"); Sys.time()

Sys.time(); trajectory_function("anisoEVI", "swindow_periodic", "central", "value_int", "value_int"); Sys.time()
Sys.time(); trajectory_function("anisoEVI", "swindow_periodic", "eastern", "value_int", "value_int"); Sys.time()
Sys.time(); trajectory_function("anisoEVI", "swindow_periodic", "southern", "value_int", "value_int"); Sys.time()
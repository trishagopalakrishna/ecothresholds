#In this script, I break up the meanfill rds files of each index for each climate zone
#into many smaller rds files containing 100 unique pixels, each pixel having 240 rows (Jan 2002 - Dec 2021)

#This script is meant to be run on the server

library(tidyverse)

chunk2 <- 240*100

break_rds_function <- function (rds_file_path, output_folder_name, climate_zone){
  rds_file <- read_rds(rds_file_path)
  
  n2 <- nrow(rds_file)
  r2  <- rep(1:ceiling(n2/chunk2),each=chunk2)[1:(n2)]
  split_pivot_df <- split(rds_file, r2)
  
  for (i in 1:length(split_pivot_df)){
    print (paste0("Counter ", i))
    write_rds(split_pivot_df[[i]], 
              paste0(file_path, output_folder_name, "split_", climate_zone, "_", i,".rds"))
  }
}

file_path <- "/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/anisoEVI/Max_consecutive_missingNA/" #manually changed index_name
file_list <- list.files(path = file_path ,pattern='.rds$', all.files=TRUE, full.names=TRUE) 

break_rds_function(file_list[[1]], "meanfill_central_dfs/", "central")
break_rds_function(file_list[[2]], "meanfill_eastern_dfs/", "eastern")
break_rds_function(file_list[[3]], "meanfill_southern_dfs/", "southern")
#In this script, I bring together the multiple rds files which are the output of 1_stl.Rmd
#into one rds file per index, per swin parameter and per climate zone

library(tidyverse)

index_file_path <- "/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/"

bind_files_function <- function (index_name, swinfolder, climate_zone){
  file_list<- list.files (path = paste0(index_file_path, index_name,"/", "STL_Decomposition", "/", swinfolder, "/"), pattern =paste0("*",climate_zone,"*"), all.files = T, full.names = T)
  file_list <- gtools::mixedsort(file_list)
  rds_list <- lapply(file_list, read_rds)
  message("completed reading all rds files")
  
  x <- bind_rows(rds_list)
  message ("completed binding to one df")
 
  write_rds(x, paste0(index_file_path, index_name, "/", "STL_Decomposition","/", swinfolder, "/", "stl_", climate_zone,".rds"))
  message ("completed writing to disk")
}


#Sys.time(); bind_files_function("NDVI", "swindow_7","central"); Sys.time()
Sys.time(); bind_files_function("NDVI", "swindow_7","eastern"); Sys.time()
Sys.time(); bind_files_function("NDVI", "swindow_7","southern"); Sys.time()

Sys.time(); bind_files_function("NDVI", "swindow_periodic","central"); Sys.time()
Sys.time(); bind_files_function("NDVI", "swindow_periodic","eastern"); Sys.time()
Sys.time(); bind_files_function("NDVI", "swindow_periodic","southern"); Sys.time()

Sys.time(); bind_files_function("EVI", "swindow_7","central"); Sys.time()
Sys.time(); bind_files_function("EVI", "swindow_7","eastern"); Sys.time()
Sys.time(); bind_files_function("EVI", "swindow_7","southern"); Sys.time()

Sys.time(); bind_files_function("EVI", "swindow_periodic","central"); Sys.time()
Sys.time(); bind_files_function("EVI", "swindow_periodic","eastern"); Sys.time()
Sys.time(); bind_files_function("EVI", "swindow_periodic","southern"); Sys.time()




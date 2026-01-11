library(tidyverse)
library(foreach)
library(doParallel)

numCores<- 16
my.cluster<- parallel::makeCluster(
  numCores,
  type = "FORK",
  outfile = ""
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered() #should be TRUE
foreach::getDoParWorkers() #should give numCores

index_file_path <- "/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/"

stl_function <- function (index_name, input_rds_climate_zone_folder, swin, swinfolder, climate_zone){
  file_list<- list.files (path = paste0(index_file_path, index_name,"/", "Max_consecutive_missingNA/", input_rds_climate_zone_folder, "/"), pattern =".rds$", all.files = T, full.names = T)
  file_list <- gtools::mixedsort(file_list)
  x<- split(file_list, ceiling(seq_along(file_list)/1000))
  
  for (i in 1:length(x)){
    rds_list<- lapply(x[[i]], read_rds)
    stl_decomposition<- function (df, swin){
      df_ts<- ts(df$value_int, start=c(2002,1), end=c(2021,12), frequency=12)
      stl_df<- stlplus::stlplus(df_ts, s.window = swin, s.degree = 1, t.degree = 1)
      df_needed<- (stl_df$data)
      df_needed<- df_needed %>% dplyr::select(-c(weights, sub.labels ))
      df_needed
    }
    message ("Stl function defined")
    
    results <- foreach(
      i= 1:length(rds_list),
      .combine = "rbind") %dopar% {
        rds_list[[i]] %>%
          group_by(cell) %>%
          nest() %>%
          mutate(stl_data = purrr::map(data, stl_decomposition, swin= swin)) %>%
          unnest(c(data,stl_data)) %>%
          dplyr::select(-c(Year, Month, raw))
      }
    message("Stl decomposition complete")
    
    write_rds(results, paste0(index_file_path, index_name, "/", "STL_Decomposition","/", swinfolder,"/", "stl_", climate_zone, "_", i, ".rds"))#changed output file address manually
    message ("rds written out")
  }
}

Sys.time(); stl_function("NDVI", "meanfill_central_dfs", 7, "swindow_7", "central"); Sys.time()
Sys.time(); stl_function("NDVI", "meanfill_eastern_dfs", 7, "swindow_7", "eastern"); Sys.time()
Sys.time(); stl_function("NDVI", "meanfill_southern_dfs", 7, "swindow_7", "southern"); Sys.time()

Sys.time(); stl_function("EVI", "meanfill_central_dfs", 7, "swindow_7", "central"); Sys.time()
Sys.time(); stl_function("EVI", "meanfill_eastern_dfs", 7, "swindow_7", "eastern"); Sys.time()
Sys.time(); stl_function("EVI", "meanfill_southern_dfs", 7, "swindow_7", "southern"); Sys.time()

Sys.time(); stl_function("NDVI", "meanfill_central_dfs", "periodic", "swindow_periodic", "central"); Sys.time()
Sys.time(); stl_function("NDVI", "meanfill_eastern_dfs", "periodic", "swindow_periodic", "eastern"); Sys.time()
Sys.time(); stl_function("NDVI", "meanfill_southern_dfs", "periodic", "swindow_periodic", "southern"); Sys.time()

Sys.time(); stl_function("EVI", "meanfill_central_dfs", "periodic", "swindow_periodic", "central"); Sys.time()
Sys.time(); stl_function("EVI", "meanfill_eastern_dfs", "periodic", "swindow_periodic", "eastern"); Sys.time()
Sys.time(); stl_function("EVI", "meanfill_southern_dfs", "periodic", "swindow_periodic", "southern"); Sys.time()


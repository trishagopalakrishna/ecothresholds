library(tidyverse)
library(foreach)
library(doParallel)


model_levels <- c("Null", "Lin", "Step", "Quad")

model_select<- function(trajectories_df){ 
  trajectories_df <- trajectories_df %>% mutate(model_order=ordered(c("Quad", "Lin", "Null", "Step"),
                                                                    levels = model_levels))
  
  message("Least AIC model selection")
  trajectories_df <- trajectories_df %>%
    mutate(aic_diff= aic - min(aic)) #wrt model with least AIC
  
  message("Conditional model selection")
  condition_less2<- trajectories_df %>% filter(aic_diff<=2)
  
  if (dim(condition_less2)[1]==1){
    trajectories_df<- condition_less2 
  } else{
    trajectories_df<- condition_less2 %>% 
      filter(model_order==min(model_order))
  }
  trajectories_df 
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

index_file_path <- "/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/"

model_selection_function <- function (index_name, swinfolder_input, annual_or_monthly, swinfolder_output) {
  file_list<- list.files (path = paste0(index_file_path, index_name,"/", "TrajectoryShapes/", annual_or_monthly, "/", swinfolder_input, "/"), 
                          pattern =".rds$", all.files = T, full.names = T)
  file_list <- gtools::mixedsort(file_list)
  message ("reading file paths complete")
  
  y<- foreach(
    i = 1:length(file_list)) %dopar% {
      read_rds(file_list[[i]]) %>%
        group_by(cell,x,y) %>%
        nest() %>%  
        mutate(selection_data = purrr::map(data, model_select)) %>%
        select(-data) %>%
        unnest(selection_data) %>%
        mutate(climate_zone = str_split(str_split(str_split(file_list[[i]], "/")[[1]][14], "_")[[1]][2], ".rds")[[1]][1])
    }
  compiled_df <- bind_rows(y)
  write_rds(compiled_df, paste0(index_file_path, index_name, "/", "TrajectoryShapes/", annual_or_monthly,"/", swinfolder_output, "/", "modelselection_", swinfolder_input, ".rds"))
  
}

#NDVI
Sys.time(); model_selection_function ("NDVI", "swindow_7", "monthly", "swindow_7"); Sys.time()
Sys.time(); model_selection_function ("NDVI", "swindow_11", "monthly", "swindow_11"); Sys.time()
Sys.time(); model_selection_function ("NDVI", "swindow_periodic","monthly", "swindow_periodic"); Sys.time()
Sys.time(); model_selection_function ("NDVI", "value_int","monthly", "value_int"); Sys.time()

#EVI
#Sys.time(); model_selection_function ("EVI", "swindow_7","monthly", "swindow_7"); Sys.time()
#Sys.time(); model_selection_function ("EVI", "swindow_11","monthly", "swindow_11"); Sys.time()
#Sys.time(); model_selection_function ("EVI", "swindow_periodic", "monthly", "swindow_periodic"); Sys.time()
#Sys.time(); model_selection_function ("EVI", "value_int", "monthly", "value_int"); Sys.time()

#anisoEVI
#Sys.time(); model_selection_function ("anisoEVI", "swindow_7", "monthly", "swindow_7"); Sys.time()
#Sys.time(); model_selection_function ("anisoEVI", "swindow_11", "monthly", "swindow_11"); Sys.time()
#Sys.time(); model_selection_function ("anisoEVI", "swindow_periodic", "monthly" ,"swindow_periodic"); Sys.time()
#Sys.time(); model_selection_function ("anisoEVI", "value_int", "monthly" ,"value_int"); Sys.time()

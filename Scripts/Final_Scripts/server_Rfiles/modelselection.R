library(tidyverse)
library(foreach)
library(doParallel)



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
  
  if (dim(condition_less2)[1]==1){
    models_pixel_df<- condition_less2 
  } else{
    models_pixel_df<- condition_less2 %>% 
      filter(model_order==min(model_order))
  }
  models_pixel_df 
}

rds_list <- list.files(path = "/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/anisoEVI/TrajectoryShapes/annual/swindow_7/", 
                       pattern='annuals7_trajshape_southern', all.files=TRUE, full.names=TRUE) #manually changed to read all files of each region

data_rds_list<- list()
for (i in 1:length(rds_list)){
  data_rds_list[[i]]<- read_rds(rds_list[i])
}

data_rds_list<- data_rds_list %>% discard(is.null)

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


Sys.time(); y<- foreach (
  i=1:length(data_rds_list),
  .combine = "rbind") %dopar% 
  (data_rds_list[[i]] %>%
     rename("trend_name"="trend") %>%
     dplyr::group_by(cell) %>%
     nest() %>%
     mutate(modelselect_data = purrr::map(data, model_select)) %>%
     select(-data) %>%
     unnest(modelselect_data)
  ); Sys.time()

write_rds(y, "/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/anisoEVI/TrajectoryShapes/annual/swindow_7/finalshape_annuals7_southern.rds") #manually changed to make one file per region

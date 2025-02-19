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

Sys.time()
file_list <- list.files(path = "/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/NDVI/Max_consecutive_missingNA/meanfill_southern_dfs/", 
                        pattern='.rds$', all.files=TRUE, full.names=TRUE)
rds_list<-lapply(file_list, read_rds)
Sys.time()
length(rds_list)

Sys.time(); print ("Unique number of cells check")
unique_cell_df<- data.frame(matrix(ncol = 2, nrow = 0))
x <- c("df_number", "number_uniquecells")
colnames(unique_cell_df) <- x
for (i in 2001:3000){
  unique_cell_df[i,1]<-i
  unique_cells<-length(unique(rds_list[[i]]$cell))
  unique_cell_df[i,2]<- unique_cells
}

write.csv(unique_cell_df, "/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/NDVI/STL_Decomposition/swindow_11/unique_cell_nrowcheck.csv")

Sys.time(); print ("STL start")
stl<- function (split_pivot_df){
  stl_decomposition<- function (df, swin){
    df_ts<- ts(df$value_int, start=c(2002,1), end=c(2021,12), frequency=12)
    stl_df<- stlplus::stlplus(df_ts, s.window = swin, s.degree = 1, t.degree = 1)
    df_needed<- (stl_df$data)
    df_needed<- df_needed %>% dplyr::select(-c(weights, sub.labels ))
    df_needed
  }
  message ("Stl function defined")
  
  results <- foreach(
    i= 2001:3000,
    .combine = "rbind") %dopar% {
      split_pivot_df[[i]] %>%
        group_by(cell) %>%
        nest() %>%
        mutate(stl_data = purrr::map(data, stl_decomposition, swin=11)) %>%
        unnest(c(data,stl_data)) %>%
        dplyr::select(-c(Year, Month, raw))
    }
  message("Stl decomposition complete")
  write_rds(results, "/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/NDVI/STL_Decomposition/swindow_11/stl_southern2001_3000.rds")
  message ("rds written out")
}

Sys.time(); stl(rds_list); Sys.time()
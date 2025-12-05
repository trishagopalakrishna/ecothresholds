library(tidyverse)

southern_rds<- read_rds("/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/NDVI/Max_consecutive_missingNA/meanfill_NDVI_southern.rds")

Sys.time()

chunk2 <- 240*100
n2 <- nrow(central_rds)
r2  <- rep(1:ceiling(n2/chunk2),each=chunk2)[1:(n2)]
split_pivot_df <- split(central_rds, r2)
remove(chunk2, n2,r2)
remove(central_rds)
gc()

Sys.time()

for (i in 1:length(split_pivot_df)){
  print (paste0("Counter ", i))
  write_rds(split_pivot_df[[i]], 
            paste0("/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/NDVI/Max_consecutive_missingNA/meanfill_southern_dfs/", "split_southern_", i,".rds"))
}

Sys.time()
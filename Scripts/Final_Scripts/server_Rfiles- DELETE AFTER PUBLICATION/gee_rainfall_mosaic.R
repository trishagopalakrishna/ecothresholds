library(terra)

## 1- Mosaic
raster_filepath <- list.files(path = "/dungbeetle/home/tg505/Trisha/ecothresholds/Data/Climate/", pattern= ".tif$", all.files=TRUE, full.names=TRUE)
print (raster_filepath)
raster_list<-lapply(raster_filepath, rast)
print ("Mosaicing starting")
Sys.time(); gee_chirps_1km <- do.call(mosaic,raster_list); Sys.time()
print ("Writing out mosiac rainfall 1km")
writeRaster(gee_chirps_1km,  "/dungbeetle/home/tg505/Trisha/ecothresholds/Data/Climate/mosiac_gee_chirps_1km.tif")

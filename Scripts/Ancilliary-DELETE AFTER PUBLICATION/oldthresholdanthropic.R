
Copy pasted below from old 3_prep_indices.Rmd script because I cannot figure out why I did this and I donot know yet
if I need to retain this script for other purposes

# #10% anthropic mask
# thresh10_anthropic<- rast(here("Data", "Masks_to_define_pixelsofinterest", "Mapbiomas_anthropic_mask", "thresholded_mapbiomas_mask_10.tif"))
# thresh20_anthropic<- rast(here("Data", "Masks_to_define_pixelsofinterest", "Mapbiomas_anthropic_mask", "thresholded_mapbiomas_mask_20.tif"))

d_trans<- st_read(here ("Data", "Cattelanetal_Clustering", "cerrado_climate_zones.shp"))
d_trans<- st_transform(d_trans, crs = 4326)
cerrado<- d_trans %>% st_union()

# Sys.time();thresh10_central<- terra::crop(thresh10_anthropic, vect(d_trans$geometry[[1]]))
# thresh10_central_mask<- terra::mask(thresh10_central, vect(d_trans$geometry[[1]])); Sys.time()
# Sys.time(); thresh10_southern<- terra::crop(thresh10_anthropic, vect(d_trans$geometry[[2]]))
# thresh10_southern_mask<- terra::mask(thresh10_southern, vect(d_trans$geometry[[2]])); Sys.time()
# Sys.time();thresh10_eastern<- terra::crop(thresh10_anthropic, vect(d_trans$geometry[[3]]))
# thresh10_eastern_mask<- terra::mask(thresh10_eastern, vect(d_trans$geometry[[3]])); Sys.time()
# remove(thresh10_central, thresh10_southern, thresh10_eastern)
# writeRaster(thresh10_central_mask, here("Data", "Masks_to_define_pixelsofinterest", "Mapbiomas_anthropic_mask", "thresh10_central.tif"))
# writeRaster(thresh10_southern_mask, here("Data", "Masks_to_define_pixelsofinterest", "Mapbiomas_anthropic_mask", "thresh10_southern.tif"))
# writeRaster(thresh10_eastern_mask, here("Data", "Masks_to_define_pixelsofinterest", "Mapbiomas_anthropic_mask", "thresh10_eastern.tif"))
# 
# Sys.time();thresh20_central<- terra::crop(thresh20_anthropic, vect(d_trans$geometry[[1]]))
# thresh20_central_mask<- terra::mask(thresh20_central, vect(d_trans$geometry[[1]])); Sys.time()
# Sys.time(); thresh20_southern<- terra::crop(thresh20_anthropic, vect(d_trans$geometry[[2]]))
# thresh20_southern_mask<- terra::mask(thresh20_southern, vect(d_trans$geometry[[2]])); Sys.time()
# Sys.time();thresh20_eastern<- terra::crop(thresh20_anthropic, vect(d_trans$geometry[[3]]))
# thresh20_eastern_mask<- terra::mask(thresh20_eastern, vect(d_trans$geometry[[3]])); Sys.time()
# remove(thresh20_central, thresh20_southern, thresh20_eastern)
# writeRaster(thresh20_central_mask, here("Data", "Masks_to_define_pixelsofinterest", "Mapbiomas_anthropic_mask", "thresh20_central.tif"))
# writeRaster(thresh20_southern_mask, here("Data", "Masks_to_define_pixelsofinterest", "Mapbiomas_anthropic_mask", "thresh20_southern.tif"))
# writeRaster(thresh20_eastern_mask, here("Data", "Masks_to_define_pixelsofinterest", "Mapbiomas_anthropic_mask", "thresh20_eastern.tif"))
# 
# remove(thresh10_anthropic, thresh10_central_mask, thresh10_eastern_mask, thresh10_southern_mask,
#        thresh20_anthropic, thresh20_central_mask, thresh20_eastern_mask, thresh20_southern_mask)
# remove(cerrado)
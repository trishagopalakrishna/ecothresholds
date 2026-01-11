library(RColorBrewer)
library(lattice)
library(tidyverse)
library(here)
library(tictoc)


library(ggplot2)
library(ggpubr)
library(sf)
library(terra)

library(tmap)
library(tmaptools)
library(RColorBrewer)
library(viridis)

terraOptions(memfrac=0.5, tempdir = here("Scratch"), progress=10)


## In this script I prepare the training data used in Lewisetal2023 about the physiognomies in Chapada NP. I prepare this data to 
## use it as field 'test' data to better understand anisoEVI and kNDVI.
## Kennedy sent me the .kmz file on Tue 30th Jan 2024. The kmz file has individual kmls of the different vegetation formations i.e 
## GPS locations of the different types of physiognomies. I extracted individual shps of the different types of physiognomies
## in QGIS (added to map and exported as shp).

##Load individual physiognomy shps
## Non physiognomy points
ag<- st_read(here("Data", "Lewisetal2022_ChapadaNP_PhysiognomyMaps", "agriculture.shp"))
ag<- ag %>% mutate(Name= "agriculture", Formation= "Agriculture")
ag<- ag %>% dplyr::select(c(Name, Formation))
water <- st_read(here("Data", "Lewisetal2022_ChapadaNP_PhysiognomyMaps", "water.shp"))
water<- water %>% mutate(Name= "water", Formation= "Water")
water<- water %>% dplyr::select(c(Name, Formation))
plantation <- st_read(here("Data", "Lewisetal2022_ChapadaNP_PhysiognomyMaps", "plantation.shp"))
plantation<- plantation %>% mutate(Name= "plantation", Formation= "Plantation")
plantation<- plantation %>% dplyr::select(c(Name, Formation))
pasture <- st_read(here("Data", "Lewisetal2022_ChapadaNP_PhysiognomyMaps", "pasture.shp"))
pasture<- pasture %>% mutate(Name= "pasture", Formation="Pasture")
pasture<- pasture %>% dplyr::select(c(Name, Formation))
nonveg <- st_read(here("Data", "Lewisetal2022_ChapadaNP_PhysiognomyMaps", "un_veg.shp"))
nonveg<- nonveg %>% mutate(Name= "nonveg", Formation="Non Veg")
nonveg<- nonveg %>% dplyr::select(c(Name, Formation))

##Veg points
openwetgrassland<- st_read(here("Data", "Lewisetal2022_ChapadaNP_PhysiognomyMaps", "openwet_grassland.shp"))
openwetgrassland<- openwetgrassland %>% mutate(Name="campo_limpo_umido", Formation="Grassland")
openwetgrassland<- openwetgrassland %>% dplyr::select(c(Name, Formation))
camp_seco<- st_read(here("Data", "Lewisetal2022_ChapadaNP_PhysiognomyMaps", "campo_seco.shp"))
camp_seco<- camp_seco %>% mutate(Name= "camp_seco", Formation="Grassland")
camp_seco<- camp_seco %>% dplyr::select(c(Name, Formation))
camp_sujo<- st_read(here("Data", "Lewisetal2022_ChapadaNP_PhysiognomyMaps", "campo_sujo.shp"))
camp_sujo<- camp_sujo %>% mutate(Name= "camp_sujo", Formation="Grassland")
camp_sujo<- camp_sujo %>% dplyr::select(c(Name, Formation))

#Kennedy's email said that for her paper she combined campo_sujo and campo_seco, which is what I do below
grassland<- dplyr::bind_rows(camp_seco, camp_sujo)
remove(camp_seco, camp_sujo)

rupestrian_grassland<- st_read(here("Data", "Lewisetal2022_ChapadaNP_PhysiognomyMaps", "ruperstrian_grassland.shp"))
rupestrian_grassland<- rupestrian_grassland %>% mutate(Name= "campo_rupestrian", Formation="Grassland")
rupestrian_grassland<- rupestrian_grassland %>% dplyr::select(c(Name, Formation))

typical_cerrado<- st_read(here("Data", "Lewisetal2022_ChapadaNP_PhysiognomyMaps", "typical_cerrado.shp"))
typical_cerrado<- typical_cerrado %>% mutate(Name="cerrado_sensu_stricto", Formation="Savanna")
typical_cerrado<- typical_cerrado%>% dplyr::select(c(Name, Formation))
rupestrian_cerrado<- st_read(here("Data", "Lewisetal2022_ChapadaNP_PhysiognomyMaps", "ruperstrian_cerrado.shp"))
rupestrian_cerrado<- rupestrian_cerrado %>% mutate(Name="cerrado_rupestre", Formation="Savanna")
rupestrian_cerrado<- rupestrian_cerrado %>% dplyr::select(c(Name, Formation))
palm_swamp<- st_read(here("Data", "Lewisetal2022_ChapadaNP_PhysiognomyMaps", "palm_swamp.shp"))
palm_swamp<- palm_swamp %>% mutate(Name="vereda", Formation="Savanna")
palm_swamp<- palm_swamp %>% dplyr::select(c(Name, Formation))

dense_cerrado_woodland<- st_read(here("Data", "Lewisetal2022_ChapadaNP_PhysiognomyMaps", "dense_cerrado_woodland.shp"))
dense_cerrado_woodland<- dense_cerrado_woodland %>% mutate(Name="cerrado_sensu_stricto", Formation="Forest and Woodland")
dense_cerrado_woodland<- dense_cerrado_woodland %>% dplyr::select(c(Name, Formation))
gallery_forest<- st_read(here("Data", "Lewisetal2022_ChapadaNP_PhysiognomyMaps", "gallery_forest.shp"))
gallery_forest<- gallery_forest %>% mutate(Name="mata_de_galeria", Formation="Forest and Woodland")
gallery_forest<- gallery_forest %>% dplyr::select(c(Name, Formation))

lewis_points<- dplyr::bind_rows(ag, water, plantation, pasture, nonveg, openwetgrassland, grassland, rupestrian_grassland, 
                                typical_cerrado, rupestrian_cerrado, palm_swamp, dense_cerrado_woodland, gallery_forest )
lewis_points<- lewis_points %>% st_zm(drop=T)
st_write(lewis_points, here("Data", "Lewisetal2022_ChapadaNP_PhysiognomyMaps", "prepped_points.shp"))

library(tidyverse)
library(ggplot2)


swin11_southern<-read_rds("/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/NDVI/TrajectoryShapes/annual/swindow_11/finalshape_annuals11_southern.rds")
swin11_eastern<- read_rds("/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/NDVI/TrajectoryShapes/annual/swindow_11/finalshape_annuals11_eastern.rds")
swin11_central<- read_rds("/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/NDVI/TrajectoryShapes/annual/swindow_11/finalshape_annuals11_central.rds")

finalmodel_df<- rbind(swin11_southern, swin11_eastern, swin11_central)
 
detrended_final_model<- finalmodel_df %>% 
   dplyr::select(c(cell,x,y, model_order, shape_class, trend_name, loc_brk)) %>%
   rename(c("cell"="cell","x"="x", "y"="y", "Detrended_ModelOrder"="model_order", "Detrended_ShapeClass"="shape_class",
            "Detrended_TrendName"="trend_name", "Detrended_LocBrk"="loc_brk"))
 detrended_final_model<- detrended_final_model %>% mutate(DetrendedResults = paste0(Detrended_ModelOrder,"_", Detrended_ShapeClass,"_", Detrended_TrendName))
 detrended_final_model2<- detrended_final_model %>% dplyr::select(c(cell,DetrendedResults, Detrended_LocBrk))
 pivot_df<- detrended_final_model2 %>% pivot_longer(2) %>% dplyr::select(-c(name))
 
 result1a_plot_df<- pivot_df %>% group_by(value) %>% 
   summarise(count= n()) %>%
   mutate(percentage= (count/nrow(finalmodel_df))*100) 
 names(result1a_plot_df)<- c("TrajShape", "STL11Count", "STL11Percentage")

 nb.cols <- 10
 mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)
 
 r1a<- result1a_plot_df %>% 
   ggplot(aes(TrajShape, value, fill = TrajShape)) +
   geom_col(position = position_dodge(), color = "black") +
   theme_classic() +
   #facet_grid(~value, scales = "free_x", switch = "x") +
   xlab("Trajectory Shape")+ ylab ("% of total area of Cerrado")+
   theme(strip.placement  = "outside",
         panel.spacing    = unit(0, "points"),
         strip.background = element_blank(),
         strip.text       = element_text(size = 12, angle = 90),
         axis.text.x = element_text(size = 12, angle = 90))+
   xlab("Model Type")
 r1a<- r1a + scale_fill_manual(values=mycolors)
 ggsave("/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/NDVI/TrajectoryShapes/annual/swindow_11/ndvionly_s11_alltrajresults_barplot.png",
        r1a,dpi = 700, height = 20, width=20, units = "cm")
library(ggplot2)
library(tidyverse)

##Data input

#1. No stl
nostl_central<- read_rds("/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/anisoEVI/TrajectoryShapes/annual/no_detrending/finalshape_annualmean_central.rds")
nostl_southern<- read_rds("/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/anisoEVI/TrajectoryShapes/annual/no_detrending/finalshape_annualmean_southern.rds")
nostl_eastern<- read_rds("/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/anisoEVI/TrajectoryShapes/annual/no_detrending/finalshape_annualmean_eastern.rds")

#2. swin=11
swin11_central<- read_rds("/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/anisoEVI/TrajectoryShapes/annual/swindow_11/finalshape_annuals11_central.rds")
swin11_southern<- read_rds("/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/anisoEVI/TrajectoryShapes/annual/swindow_11/finalshape_annuals11_southern.rds")
swin11_eastern<- read_rds("/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/anisoEVI/TrajectoryShapes/annual/swindow_11/finalshape_annuals11_eastern.rds")

#3. swin=7
swin7_central<- read_rds("/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/anisoEVI/TrajectoryShapes/annual/swindow_7/finalshape_annuals7_central.rds")
swin7_southern<- read_rds("/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/anisoEVI/TrajectoryShapes/annual/swindow_7/finalshape_annuals7_southern.rds")
swin7_eastern<- read_rds("/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/anisoEVI/TrajectoryShapes/annual/swindow_7/finalshape_annuals7_eastern.rds")

#4. swin=periodic
swinperiodic_central<- read_rds("/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/anisoEVI/TrajectoryShapes/annual/swindow_periodic/finalshape_annualsperiodic_central.rds")
swinperiodic_southern<- read_rds("/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/anisoEVI/TrajectoryShapes/annual/swindow_periodic/finalshape_annualsperiodic_southern.rds")
swinperiodic_eastern<- read_rds("/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/anisoEVI/TrajectoryShapes/annual/swindow_periodic/finalshape_annualsperiodic_eastern.rds")

swin11_df<- rbind(swin11_southern, swin11_eastern, swin11_central)
swin7_df<- rbind(swin7_southern, swin7_eastern, swin7_central)
swinperiodic_df<- rbind(swinperiodic_southern, swinperiodic_eastern, swinperiodic_central)
nostl_df<- rbind(nostl_southern, nostl_eastern, nostl_central)

remove(nostl_central, nostl_southern, nostl_eastern, swin11_central, swin11_southern, swin11_eastern, swin7_central, swin7_southern, swin7_eastern, swinperiodic_central, swinperiodic_southern, swinperiodic_eastern)

print ("Data read")

## Data processing

#1. no stl
nostl_final_model<- nostl_df %>% 
  dplyr::select(c(cell,x,y, model_order, shape_class, trend_name, loc_brk)) %>%
  rename(c("cell"="cell","x"="x", "y"="y", "NoSTL_ModelOrder"="model_order", "NoSTL_ShapeClass"="shape_class",
           "NoSTL_TrendName"="trend_name", "NoSTL_LocBrk"="loc_brk"))
nostl_final_model<- nostl_final_model %>% mutate(NoSTLResults = paste0(NoSTL_ModelOrder,"_", NoSTL_ShapeClass,"_", NoSTL_TrendName))
nostl_final_model2<- nostl_final_model %>% dplyr::select(c(cell, NoSTLResults, NoSTL_LocBrk))
nostl_pivot_df<- nostl_final_model2 %>% pivot_longer(2) %>% dplyr::select(-c(name))

#2. swin=11
swin11_final_model<- swin11_df %>% 
  dplyr::select(c(cell,x,y, model_order, shape_class, trend_name, loc_brk)) %>%
  rename(c("cell"="cell","x"="x", "y"="y", "Swin11_ModelOrder"="model_order", "Swin11_ShapeClass"="shape_class",
           "Swin11_TrendName"="trend_name", "Swin11_LocBrk"="loc_brk"))
swin11_final_model<- swin11_final_model %>% mutate(Swin11Results = paste0(Swin11_ModelOrder,"_", Swin11_ShapeClass,"_", Swin11_TrendName))
swin11_final_model2<- swin11_final_model %>% dplyr::select(c(cell,Swin11Results, Swin11_LocBrk))
pivot_swin11_df<- swin11_final_model2 %>% pivot_longer(2) %>% dplyr::select(-c(name))

#3.swin=7
swin7_final_model<- swin7_df %>% 
  dplyr::select(c(cell,x,y, model_order, shape_class, trend_name, loc_brk)) %>%
  rename(c("cell"="cell","x"="x", "y"="y", "Swin7_ModelOrder"="model_order", "Swin7_ShapeClass"="shape_class",
           "Swin7_TrendName"="trend_name", "Swin7_LocBrk"="loc_brk"))
swin7_final_model<- swin7_final_model %>% mutate(Swin7Results = paste0(Swin7_ModelOrder,"_", Swin7_ShapeClass,"_", Swin7_TrendName))
swin7_final_model2<- swin7_final_model %>% dplyr::select(c(cell,Swin7Results, Swin7_LocBrk))
pivot_swin7_df<- swin7_final_model2 %>% pivot_longer(2) %>% dplyr::select(-c(name))

#4. swin=periodic
swinperiodic_final_model<- swinperiodic_df %>% 
  dplyr::select(c(cell,x,y, model_order, shape_class, trend_name, loc_brk)) %>%
  rename(c("cell"="cell","x"="x", "y"="y", "Swinperiodic_ModelOrder"="model_order", "Swinperiodic_ShapeClass"="shape_class",
           "Swinperiodic_TrendName"="trend_name", "Swinperiodic_LocBrk"="loc_brk"))
swinperiodic_final_model<- swinperiodic_final_model %>% mutate(SwinperiodicResults = paste0(Swinperiodic_ModelOrder,"_", Swinperiodic_ShapeClass,"_", Swinperiodic_TrendName))
swinperiodic_final_model2<- swinperiodic_final_model %>% dplyr::select(c(cell,SwinperiodicResults, Swinperiodic_LocBrk))
pivot_swinperiodic_df<- swinperiodic_final_model2 %>% pivot_longer(2) %>% dplyr::select(-c(name))

print ("Data prep 1 complete")

#barplot data prep 
nostl_plot_df<- nostl_pivot_df %>% group_by(value) %>% 
  summarise(count= n()) %>%
  mutate(percentage= (count/nrow(nostl_df))*100) 
names(nostl_plot_df)<- c("TrajShape", "AM_Count", "AM_Percentage")

swin11_plot_df<- pivot_swin11_df %>% group_by(value) %>% 
  summarise(count= n()) %>%
  mutate(percentage= (count/nrow(swin11_df))*100) 
names(swin11_plot_df)<- c("TrajShape", "AMSwin11_Count", "AMSwin11_Percentage")

swin7_plot_df<- pivot_swin7_df %>% group_by(value) %>% 
  summarise(count= n()) %>%
  mutate(percentage= (count/nrow(swin7_df))*100) 
names(swin7_plot_df)<- c("TrajShape", "AMSwin7_Count", "AMSwin7_Percentage")

swinperiodic_plot_df<- pivot_swinperiodic_df %>% group_by(value) %>% 
  summarise(count= n()) %>%
  mutate(percentage= (count/nrow(swinperiodic_df))*100) 
names(swinperiodic_plot_df)<- c("TrajShape", "AMSwinperiodic_Count", "AMSwinperiodic_Percentage")

am_allclasses_df<- full_join(nostl_plot_df, swin11_plot_df, by=join_by(TrajShape))
am_allclasses_df<- full_join(am_allclasses_df, swin7_plot_df, by=join_by(TrajShape))
am_allclasses_df<- full_join(am_allclasses_df, swinperiodic_plot_df, by=join_by(TrajShape))
pivot_am_allclasses_df<- am_allclasses_df %>% dplyr::select(c("TrajShape", "AM_Percentage", "AMSwin11_Percentage", "AMSwin7_Percentage", "AMSwinperiodic_Percentage")) %>%
  pivot_longer(2:5)

Set1palette<- brewer.pal(n=8, "Set1")

r3a<- pivot_am_allclasses_df %>% 
  ggplot(aes(TrajShape, value, fill = TrajShape, group=name, alpha=as.factor(name))) +
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
r3a<- r3a + scale_fill_manual(values=Set1palette)
ggsave("/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/anisoEVI/TrajectoryShapes/annual/comparison_plot.png",
       r3a,dpi = 700, height = 20, width=20, units = "cm")

remove(r3a, nostl_plot_df, swin11_plot_df, swin7_plot_df, swinperiodic_plot_df, am_allclasses_df, Set1palette)
#################################### Aggregating no trend classes
nostl_plot_df<- nostl_pivot_df %>% group_by(value) %>% 
  summarise(count= n()) 
nostl_plot_df<- nostl_plot_df %>% mutate(value=case_when(value=="Lin_decrease_constant_NA"~"Linear decrease",
                                                         value=="Lin_increase_constant_NA"~"Linear increase",
                                                         value=="Step_NA_decrease"~"Step decrease",
                                                         value=="Step_NA_increase"~"Step increase",
                                                         TRUE~"No trend"))
nostl_plot_df<- nostl_plot_df %>% group_by(value) %>% 
  summarise(totalcount= sum(count)) %>% 
  mutate(percentage= (totalcount/nrow(nostl_df))*100)
names(nostl_plot_df)<- c("TrajShape", "NoSTL_Count", "NoSTL_Percentage")

swin11_plot_df<- pivot_swin11_df %>% group_by(value) %>% 
  summarise(count= n())
swin11_plot_df<- swin11_plot_df %>% mutate(value=case_when(value=="Lin_decrease_constant_NA"~"Linear decrease",
                                                         value=="Lin_increase_constant_NA"~"Linear increase",
                                                         value=="Step_NA_decrease"~"Step decrease",
                                                         value=="Step_NA_increase"~"Step increase",
                                                         TRUE~"No trend"))
swin11_plot_df<- swin11_plot_df %>% group_by(value) %>% 
  summarise(totalcount= sum(count)) %>% 
  mutate(percentage= (totalcount/nrow(swin11_df))*100)
names(swin11_plot_df)<- c("TrajShape", "Swin11_Count", "Swin11_Percentage")

swin7_plot_df<- pivot_swin7_df %>% group_by(value) %>% 
  summarise(count= n())
swin7_plot_df<- swin7_plot_df %>% mutate(value=case_when(value=="Lin_decrease_constant_NA"~"Linear decrease",
                                                           value=="Lin_increase_constant_NA"~"Linear increase",
                                                           value=="Step_NA_decrease"~"Step decrease",
                                                           value=="Step_NA_increase"~"Step increase",
                                                           TRUE~"No trend"))
swin7_plot_df<- swin7_plot_df %>% group_by(value) %>% 
  summarise(totalcount= sum(count)) %>% 
  mutate(percentage= (totalcount/nrow(swin7_df))*100)
names(swin7_plot_df)<- c("TrajShape", "Swin7_Count", "Swin7_Percentage")

swinperiodic_plot_df<- pivot_swinperiodic_df %>% group_by(value) %>% 
  summarise(count= n())
swinperiodic_plot_df<- swinperiodic_plot_df%>% mutate(value=case_when(value=="Lin_decrease_constant_NA"~"Linear decrease",
                                                         value=="Lin_increase_constant_NA"~"Linear increase",
                                                         value=="Step_NA_decrease"~"Step decrease",
                                                         value=="Step_NA_increase"~"Step increase",
                                                         TRUE~"No trend"))
swinperiodic_plot_df<- swinperiodic_plot_df %>% group_by(value) %>% 
  summarise(totalcount= sum(count)) %>% 
  mutate(percentage= (totalcount/nrow(swinperiodic_df))*100)
names(swinperiodic_plot_df)<- c("TrajShape", "Swinperiodic_Count", "Swinperiodic_Percentage")

am_allclasses_df<- full_join(nostl_plot_df, swin11_plot_df, by=join_by(TrajShape))
am_allclasses_df<- full_join(am_allclasses_df, swin7_plot_df, by=join_by(TrajShape))
am_allclasses_df<- full_join(am_allclasses_df, swinperiodic_plot_df, by=join_by(TrajShape)) 
pivot_am_allclasses_df<- am_allclasses_df %>% dplyr::select(c("TrajShape", "NoSTL_Percentage", "Swin11_Percentage", "Swin7_Percentage", "Swinperiodic_Percentage")) %>%
  pivot_longer(2:5)

Set1palette<- brewer.pal(n=5, "Set1")
Set1palette[[3]]<- "grey"

r3a<- pivot_am_allclasses_df %>% 
  ggplot(aes(TrajShape, value, fill = TrajShape, group=name, alpha=as.factor(name))) +
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
r3a<- r3a + scale_fill_manual(values=Set1palette)
ggsave("/dungbeetle/home/tg505/Trisha/ecothresholds/Outputs/Indices/anisoEVI/TrajectoryShapes/annual/aggtrajshapes_comparison_plot.png",
       r3a,dpi = 700, height = 20, width=20, units = "cm")



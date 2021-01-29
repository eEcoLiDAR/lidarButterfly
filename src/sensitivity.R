library(ggplot2)
library(GGally)

library(dplyr)
library(stringr)

library(gridExtra)
library(ggpubr)

setwd("D:/Sync/_Amsterdam/_PhD/")

data1=read.csv("Butterfly_lidarmetrics_25m_v3.csv")
data2=read.csv("Butterfly_lidarmetrics_25m_v3_GroningenReq.csv")

data=read.csv("Butterfly_lidarmetrics_25m_sens005.csv")

data_orig=rbind(data1,data2)

# select attributes

data_orig_sel=subset(data_orig, select=c("Transect","Transect_sec","nofret_pheightlay_b02","nofret_pheightlay_02_1","nofret_pheightlay_1_5","nofret_pheightlay_5_20","nofret_pheightlay_a20",
                                         "z090quantile.90.","rough_dsm_med","rough_dsm_med_b1"))

data_sel=subset(data, select=c("Transect","Transect_sec","nofret_pheightlay_b02","nofret_pheightlay_02_1","nofret_pheightlay_1_5","nofret_pheightlay_5_20","nofret_pheightlay_a20",
                               "z090quantile.90.","rough_dsm_med","rough_dsm_med_b1"))

data_forsens=merge(x=data_orig_sel,y=data_sel, by=c("Transect","Transect_sec"))

# visualization

p1=ggplot(data=data_forsens, aes(x=nofret_pheightlay_b02.x , y=nofret_pheightlay_b02.y),show.legend = FALSE) +  
  geom_point(color="blue",size=5,show.legend = TRUE) +
  geom_abline()+
  theme_bw(base_size = 30)+
  xlab("Density < 0.2 m [%]")+
  ylab("Modified density < 0.2 m [%]")+
  stat_cor(method = "spearman", label.x = 0, label.y = 35,size=10)

p2=ggplot(data=data_forsens, aes(x=nofret_pheightlay_02_1.x , y=nofret_pheightlay_02_1.y),show.legend = FALSE) +  
  geom_point(color="blue",size=5,show.legend = TRUE) +
  geom_abline()+
  theme_bw(base_size = 30)+
  xlab("Density < 1 m [%]")+
  ylab("Modified density < 1 m [%]")+
  stat_cor(method = "spearman", label.x = 0, label.y = 0,size=10)

p3=ggplot(data=data_forsens, aes(x=nofret_pheightlay_1_5.x , y=nofret_pheightlay_1_5.y),show.legend = FALSE) +  
  geom_point(color="blue",size=5,show.legend = TRUE) +
  geom_abline()+
  theme_bw(base_size = 30)+
  xlab("Density < 5 m [%]")+
  ylab("Modified density < 5 m [%]")+
  stat_cor(method = "spearman", label.x = 0, label.y = 0,size=10)

p4=ggplot(data=data_forsens, aes(x=nofret_pheightlay_5_20.x , y=nofret_pheightlay_5_20.y),show.legend = FALSE) +  
  geom_point(color="blue",size=5,show.legend = TRUE) +
  geom_abline()+
  theme_bw(base_size = 30)+
  xlab("Density < 20 m [%]")+
  ylab("Modified density < 20 m [%]")+
  stat_cor(method = "spearman", label.x = 0, label.y = 0,size=10)

p5=ggplot(data=data_forsens, aes(x=nofret_pheightlay_a20.x , y=nofret_pheightlay_a20.y),show.legend = FALSE) +  
  geom_point(color="blue",size=5,show.legend = TRUE) +
  geom_abline()+
  theme_bw(base_size = 30)+
  xlab("Density > 20 m [%]")+
  ylab("Modified density > 20 m [%]")+
  stat_cor(method = "spearman", label.x = 0, label.y = 50,size=10)

p6=ggplot(data=data_forsens, aes(x=z090quantile.90..x , y=z090quantile.90..y),show.legend = FALSE) +  
  geom_point(color="blue",size=5,show.legend = TRUE) +
  geom_abline()+
  theme_bw(base_size = 30)+
  xlab("Height [m]")+
  ylab("Modified height [m]")+
  stat_cor(method = "spearman", label.x = 0, label.y = 33,size=10)+
  xlim(0,35)+ylim(0,35)

p7=ggplot(data=data_forsens, aes(x=rough_dsm_med.x , y=rough_dsm_med.y),show.legend = FALSE) +  
  geom_point(color="blue",size=5,show.legend = TRUE) +
  geom_abline()+
  theme_bw(base_size = 30)+
  xlab("Total veg. roughness [m]")+
  ylab("Modified total veg. roughness [m]")+
  stat_cor(method = "spearman", label.x = 0, label.y = 13,size=10)

p8=ggplot(data=data_forsens, aes(x=rough_dsm_med_b1.y , y=rough_dsm_med_b1.x),show.legend = FALSE) +  
  geom_point(color="blue",size=5,show.legend = TRUE) +
  geom_abline()+
  theme_bw(base_size = 30)+
  xlab("Low veg. roughness [m]")+
  ylab("Modified Low veg. roughness [m]")+
  stat_cor(method = "spearman", label.x = 0, label.y = 0.45,size=10)

fig=grid.arrange(p1,p5,p6,p7,p8,
                  ncol=4,
                  nrow=3)

ggsave("Fig_sens.png",plot = fig,width = 25, height = 25)

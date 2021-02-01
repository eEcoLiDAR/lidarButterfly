library(ggplot2)
library(GGally)

library(dplyr)
library(stringr)

library(gridExtra)
library(ggpubr)

setwd("D:/Sync/_Amsterdam/_PhD/")

data1=read.csv("Boloria_metrics_0420.csv")
data2=read.csv("Hipparchia_metrics_0420.csv")

data3=read.csv("Limenitis_metrics_0420.csv")
data4=read.csv("Melitaea_metrics_0420.csv")

data25m=read.csv("Butterfly_lidarmetrics_25m_sens005.csv")
data100m=read.csv("Butterfly_lidarmetrics_100m_sens005_round2.csv")

data_orig=rbind(data1,data2,data3,data4)

# select attributes

data_orig_sel=subset(data_orig, select=c("Transect","Tr_sec","dens.b02","dens.02_1","dens.1_5","dens.5_20","dens.a20",
                                         "height","rough_total","rough_low","slope","area_low","patches_low","edges_low"))

data_sel=subset(data25m, select=c("Transect","Transect_sec","nofret_pheightlay_b02","nofret_pheightlay_02_1","nofret_pheightlay_1_5","nofret_pheightlay_5_20","nofret_pheightlay_a20",
                               "z090quantile.90.","rough_dsm_mean","rough_dsm_mean_b1"))
data_sel2=subset(data100m, select=c("Transect","Transect_sec","dtm_slope_mean","class_area_low","nofpatches_low","total_edge_low"))

data_new=merge(x=data_sel,y=data_sel2, by=c("Transect","Transect_sec"))
names(data_new)<-c("Transect","Tr_sec","dens.b02","dens.02_1","dens.1_5","dens.5_20","dens.a20",
                   "height","rough_total","rough_low","slope","area_low","patches_low","edges_low")

data_forsens=merge(x=data_orig_sel,y=data_new, by=c("Transect","Tr_sec"))

# visualization

p1=ggplot(data=data_forsens, aes(x=dens.b02.x , y=dens.b02.y),show.legend = FALSE) +  
  geom_point(color="blue",size=5,show.legend = TRUE) +
  geom_abline()+
  theme_bw(base_size = 30)+
  xlab("Density < 0.2 m [%]")+
  ylab("Modified density < 0.2 m [%]")+
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",size=10,p.accuracy = 0.001,cor.coef.name="r")

p2=ggplot(data=data_forsens, aes(x=dens.02_1.x , y=dens.02_1.y),show.legend = FALSE) +  
  geom_point(color="blue",size=5,show.legend = TRUE) +
  geom_abline()+
  theme_bw(base_size = 30)+
  xlab("Density < 1 m [%]")+
  ylab("Modified density < 1 m [%]")+
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",size=10,p.accuracy = 0.001,cor.coef.name="r")

p3=ggplot(data=data_forsens, aes(x=dens.1_5.x , y=dens.1_5.y),show.legend = FALSE) +  
  geom_point(color="blue",size=5,show.legend = TRUE) +
  geom_abline()+
  theme_bw(base_size = 30)+
  xlab("Density < 5 m [%]")+
  ylab("Modified density < 5 m [%]")+
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",size=10,p.accuracy = 0.001,cor.coef.name="r")

p4=ggplot(data=data_forsens, aes(x=dens.5_20.x , y=dens.5_20.y),show.legend = FALSE) +  
  geom_point(color="blue",size=5,show.legend = TRUE) +
  geom_abline()+
  theme_bw(base_size = 30)+
  xlab("Density < 20 m [%]")+
  ylab("Modified density < 20 m [%]")+
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",size=10,p.accuracy = 0.001,cor.coef.name="r")

p5=ggplot(data=data_forsens, aes(x=dens.a20.x , y=dens.a20.y),show.legend = FALSE) +  
  geom_point(color="blue",size=5,show.legend = TRUE) +
  geom_abline()+
  theme_bw(base_size = 30)+
  xlab("Density > 20 m [%]")+
  ylab("Modified density > 20 m [%]")+
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",size=10,p.accuracy = 0.001,cor.coef.name="r")

p6=ggplot(data=data_forsens, aes(x=height.x , y=height.y),show.legend = FALSE) +  
  geom_point(color="blue",size=5,show.legend = TRUE) +
  geom_abline()+
  theme_bw(base_size = 30)+
  xlab("Height [m]")+
  ylab("Modified height [m]")+
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",size=10,p.accuracy = 0.001,cor.coef.name="r")

p7=ggplot(data=data_forsens, aes(x=rough_total.x , y=rough_total.y),show.legend = FALSE) +  
  geom_point(color="blue",size=5,show.legend = TRUE) +
  geom_abline()+
  theme_bw(base_size = 30)+
  xlab("Total veg. roughness [m]")+
  ylab("Modified total veg. roughness [m]")+
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",size=10,p.accuracy = 0.001,cor.coef.name="r")

p8=ggplot(data=data_forsens, aes(x=rough_low.y , y=rough_low.x),show.legend = FALSE) +  
  geom_point(color="blue",size=5,show.legend = TRUE) +
  geom_abline()+
  theme_bw(base_size = 30)+
  xlab("Low veg. roughness [m]")+
  ylab("Modified low veg. roughness [m]")+
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",size=10,p.accuracy = 0.001,cor.coef.name="r")

p9=ggplot(data=data_forsens, aes(x=slope.y , y=slope.x),show.legend = FALSE) +  
  geom_point(color="blue",size=5,show.legend = TRUE) +
  geom_abline()+
  theme_bw(base_size = 30)+
  xlab("Slope [degree]")+
  ylab("Modified slope [degree]")+
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",size=10,p.accuracy = 0.001,cor.coef.name="r")

p10=ggplot(data=data_forsens, aes(x=area_low.x , y=area_low.x),show.legend = FALSE) +  
  geom_point(color="blue",size=5,show.legend = TRUE) +
  geom_abline()+
  theme_bw(base_size = 30)+
  xlab("Open area [ha]")+
  ylab("Modified open area [ha]")+
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",size=10,p.accuracy = 0.001,cor.coef.name="r")

p11=ggplot(data=data_forsens, aes(x=patches_low.x , y=patches_low.x),show.legend = FALSE) +  
  geom_point(color="blue",size=5,show.legend = TRUE) +
  geom_abline()+
  theme_bw(base_size = 30)+
  xlab("Open patches [count]")+
  ylab("Modified open patches [count]")+
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",size=10,p.accuracy = 0.001,cor.coef.name="r")

p12=ggplot(data=data_forsens, aes(x=edges_low.x , y=edges_low.x),show.legend = FALSE) +  
  geom_point(color="blue",size=5,show.legend = TRUE) +
  geom_abline()+
  theme_bw(base_size = 30)+
  xlab("Edge extent [m]")+
  ylab("Modified edge extent [m]")+
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",size=10,p.accuracy = 0.001,cor.coef.name="r")

fig=grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,
                  ncol=3,
                  nrow=4)

ggsave("Fig_sens.png",plot = fig,width = 25, height = 25)

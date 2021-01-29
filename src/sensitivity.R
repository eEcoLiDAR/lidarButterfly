library(ggplot2)
library(GGally)

library(dplyr)
library(stringr)

setwd("D:/Sync/_Amsterdam/_PhD/")

data1=read.csv("Butterfly_lidarmetrics_25m_v3.csv")
data2=read.csv("Butterfly_lidarmetrics_25m_v3_GroningenReq.csv")

data=read.csv("Butterfly_lidarmetrics_25m_sens005.csv")

data_orig=rbind(data1,data2)

# select attributes

data_orig_sel=subset(data_orig, select=c("Transect","Transect_sec","nofret_pheightlay_b02","nofret_pheightlay_02_1","nofret_pheightlay_1_5","nofret_pheightlay_5_20","nofret_pheightlay_a20",
                                         "z090quantile.90.","rough_dsm_med","rough_dsm_med_b1","dtm_slope_mean",))

data_sel=subset(data, select=c("Transect","Transect_sec","nofret_pheightlay_b02","zmean","z090quantile.90.","var_nofret",
                                         "pulsepen","zkurto","zkurto_b1","zsd","z025quantile.25.","z050quantile.50.","z075quantile.75.",
                                         "shannon"))

names(data_orig_sel) <- c("Transect","Transect_sec","Below 0.2 m density","Height","Slope")
names(data_sel) <- c("Transect","Transect_sec","Below 0.2 m density mod","Height mod","Slope mod")

data_forsens=merge(x=data_orig_selsel,y=data_selsel, by=c("Transect","Transect_sec"))

# visualization
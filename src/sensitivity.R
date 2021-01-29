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

data_orig_sel=subset(data_orig, select=c("Transect","Transect_sec","nofret_pheightlay_b02","zmean","z090quantile.90.","var_nofret",
                                         "pulsepen","zkurto","zkurto_b1","zsd","z025quantile.25.","z050quantile.50.","z075quantile.75.",
                                         "shannon"))

data_sel=subset(data, select=c("Transect","Transect_sec","nofret_pheightlay_b02","zmean","z090quantile.90.","var_nofret",
                                         "pulsepen","zkurto","zkurto_b1","zsd","z025quantile.25.","z050quantile.50.","z075quantile.75.",
                                         "shannon"))

data_orig_selsel=subset(data_orig, select=c("Transect","Transect_sec","nofret_pheightlay_b02","z090quantile.90.","dtm_slope_mean"))

data_selsel=subset(data, select=c("Transect","Transect_sec","nofret_pheightlay_b02","z090quantile.90.","dtm_slope_mean"))

names(data_orig_selsel) <- c("Transect","Transect_sec","Below 0.2 m density","Height","Slope")
names(data_selsel) <- c("Transect","Transect_sec","Below 0.2 m density mod","Height mod","Slope mod")

data_forsens=merge(x=data_orig_selsel,y=data_selsel, by=c("Transect","Transect_sec"))

#corr

tiff("corBplot.tiff", width = 6, height = 5, units = 'in', res = 600)
ggcorr(data_forsens[,c(3:8)], c("pairwise", "spearman"), name = expression(italic("Spearman's r")), label=TRUE, label_alpha=TRUE, label_size=2.5, hjust=0.85, size=3, layout.exp=3)
dev.off()
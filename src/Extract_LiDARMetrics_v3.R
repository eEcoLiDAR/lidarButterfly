"
@author: Zsofia Koma, UvA
Aim: LiDAR metrics extraction for the butterfly project version 3 with lanscape metrics
"
library(dplyr)
library(lidR)
library(rgdal)
library(raster)
library(sp)
library(sf)
library(dplyr)

library(e1071)
library(landscapemetrics)

# Set working dirctory
#workingdirectory="D:/Sync/_Amsterdam/08_coauthor_MScProjects/Reinier/lidarmetrics_calc/"
workingdirectory="D:/Reinier/Reinier_output/"

setwd(workingdirectory)
radii=100

butterflyspfile="transects.shp"

#Import
butterflysp = readOGR(dsn=butterflyspfile)
butterflysp_df=butterflysp@data

butterflysp_df_gr <- butterflysp_df %>%
  group_by(Transect) %>%
  summarise(nofobs = length(Transect))

# Direct point cloud based metrics
Transect <- as.vector(butterflysp_df_gr$Transect)
Transect= as.numeric(Transect)
#Transect=c(123)

dpcloudfea_exp_df <- data.frame(matrix(ncol = 13, nrow = 0))
x <- c("Transect", "Transect_sec","propoflowveg","class_area_low","class_area_med","class_area_high","total_edge_low","total_edge_med","total_edge_high","nofpatches_low","nofpatches_med","nofpatches_high",
       "cohesion_low","cohesion_med","cohesion_high","nofret_pheightlay_02_1","nofret_pheightlay_1_5","nofret_pheightlay_5_20")

colnames(dpcloudfea_exp_df) <- x

for (i in Transect) {
  print(i)
  
  if (file.exists(paste("Transect_",i,".laz",sep=""))) {
    
    las=readLAS(paste("Transect_",i,".laz",sep=""))
    
    butterflysp_df_sel=butterflysp_df[butterflysp_df$Transect==i,]
    
    for (j in seq(from=1,to=length(butterflysp_df_sel$x))) {
      
      las_clip=lasclipCircle(las,butterflysp_df_sel$x[j],butterflysp_df_sel$y[j],radii)
      
      if ((nrow(las_clip@data[las_clip@data$Classification==1L])>0) & (nrow(las_clip@data[las_clip@data$Classification==2L])>0)) {
        
        #writeLAS(las_clip,paste("Transect",i,"Tr_sec",butterflysp_df_sel$Tr_sec[j],"_clip_25.laz",sep=""))
        
        las_norm=lasnormalize(las_clip, knnidw(k=10,p=2))
        las_norm_veg=lasfilter(las_norm,Classification==1L)
        
        dsm = grid_metrics(las_norm,~quantile(Z, 0.90),res=1)
        proj4string(dsm) <- CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs")
        
        butterflysp_df_sel_sel=butterflysp_df_sel[butterflysp_df_sel$Tr_sec==butterflysp_df_sel$Tr_sec[j],]
        
        coordinates(butterflysp_df_sel_sel)=~x+y
        proj4string(butterflysp_df_sel_sel)<- CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs")
        
        # corrected height layers
        
        nofret_pheightlay_02_1=(nrow(las_norm@data[(las_norm@data$Classification==1L & las_norm@data$Z>0.2 & las_norm@data$Z<1)])/length(las_norm@data$Z))*100
        nofret_pheightlay_1_5=(nrow(las_norm@data[(las_norm@data$Classification==1L & las_norm@data$Z>1 & las_norm@data$Z<5)])/length(las_norm@data$Z))*100
        nofret_pheightlay_5_20=(nrow(las_norm@data[(las_norm@data$Classification==1L & las_norm@data$Z>5 & las_norm@data$Z<20)])/length(las_norm@data$Z))*100
        
        
        # landscape metrics
        
        height_class=reclassify(dsm, c(-Inf,1,1,1,5,2,5,Inf,3))
        
        class_area = sample_lsm(height_class, butterflysp_df_sel_sel,size=radii,level = "class", metric = "ca",count_boundary = FALSE,directions = 8)
        class_area_low=class_area$value[1]
        class_area_med=class_area$value[2]
        class_area_high=class_area$value[3]
        
        total_edge = sample_lsm(height_class, butterflysp_df_sel_sel,size=radii,level = "class", metric = "te",count_boundary = FALSE,directions = 8)
        total_edge_low=total_edge$value[1]
        total_edge_med=total_edge$value[2]
        total_edge_high=total_edge$value[3]
        
        nofpatches = sample_lsm(height_class, butterflysp_df_sel_sel,size=radii,level = "class", metric = "np",count_boundary = FALSE,directions = 8)
        nofpatches_low=nofpatches$value[1]
        nofpatches_med=nofpatches$value[2]
        nofpatches_high=nofpatches$value[3]
        
        cohesion = sample_lsm(height_class, butterflysp_df_sel_sel,size=radii,level = "class", metric = "cohesion",count_boundary = FALSE,directions = 8)
        cohesion_low=cohesion$value[1]
        cohesion_med=cohesion$value[2]
        cohesion_high=cohesion$value[3]
        
        # proportion of low vegetation (below 1 m)
        lowveg_class=reclassify(dsm, c(-Inf,1,1,1,Inf,NA))
        
        propoflowveg <- extract(lowveg_class,butterflysp_df_sel_sel,buffer = radii,fun=length)
        
        if (is.na(propoflowveg)) {propoflowveg=0} else {propoflowveg=propoflowveg}
        
        # Export
        
        newline <- data.frame(t(c(Transect=i,Transect_sec=paste(butterflysp_df_sel$Tr_sec[j],sep=""),
                                  propoflowveg=propoflowveg,
                                  class_area_low=class_area_low,
                                  class_area_med=class_area_med,
                                  class_area_high=class_area_high,
                                  total_edge_low=total_edge_low,
                                  total_edge_med=total_edge_med,
                                  total_edge_high=total_edge_high,
                                  nofpatches_low=nofpatches_low,
                                  nofpatches_med=nofpatches_med,
                                  nofpatches_high=nofpatches_high,
                                  cohesion_low=cohesion_low,
                                  cohesion_med=cohesion_med,
                                  cohesion_high=cohesion_high,
                                  nofret_pheightlay_02_1=nofret_pheightlay_02_1,
                                  nofret_pheightlay_1_5=nofret_pheightlay_1_5,
                                  nofret_pheightlay_5_20=nofret_pheightlay_5_20)))
        
        dpcloudfea_exp_df <- rbind(dpcloudfea_exp_df, newline)
        
      }
    }
  }
}

write.csv(dpcloudfea_exp_df,"Butterfly_lidarmetrics_100m_20200207.csv")
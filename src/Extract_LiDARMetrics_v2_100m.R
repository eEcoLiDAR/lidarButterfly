"
@author: Zsofia Koma, UvA
Aim: LiDAR metrics extraction for the butterfly project version 2
"
library(dplyr)
library(lidR)
library(rgdal)
library(raster)
library(sp)

library(e1071)

# Set working dirctory
#workingdirectory="D:/Sync/_Amsterdam/08_coauthor_MScProjects/Reinier/lidarmetrics_calc/"
workingdirectory="D:/Reinier/Reinier_output/"

setwd(workingdirectory)

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
#Transect=c(1437)

dpcloudfea_exp_df <- data.frame(matrix(ncol = 38, nrow = 0))
x <- c("Transect", "Transect_sec","rough_dsm_med","rough_dsm_med_b1","var_dsm_med","var_dsm_med_b1","var_dim_med",
       "var_dim_med_b1","dtm_slope_med","dtm_asp_med","rough_dsm_mean","rough_dsm_mean_b1","var_dsm_mean","var_dsm_mean_b1","var_dim_mean","var_dim_mean_b1","dtm_slope_mean","dtm_asp_mean","rough_dsm_sd","rough_dsm_sd_b1",
       "var_dsm_sd","var_dsm_sd_b1","var_dim_sd","var_dim_med_b1","dtm_slope_sd","dtm_asp_sd","propofhighveg")

colnames(dpcloudfea_exp_df) <- x

for (i in Transect) {
  print(i)
  
  if (file.exists(paste("Transect_",i,".laz",sep=""))) {
    
    las=readLAS(paste("Transect_",i,".laz",sep=""))
    
    butterflysp_df_sel=butterflysp_df[butterflysp_df$Transect==i,]
    
    for (j in seq(from=1,to=length(butterflysp_df_sel$x))) {
      
      las_clip=lasclipCircle(las,butterflysp_df_sel$x[j],butterflysp_df_sel$y[j],100)
      
      if ((nrow(las_clip@data[las_clip@data$Classification==1L])>0) & (nrow(las_clip@data[las_clip@data$Classification==2L])>0)) {
        
        #writeLAS(las_clip,paste("Transect",i,"Tr_sec",butterflysp_df_sel$Tr_sec[j],"_clip_25.laz",sep=""))
        
        las_norm=lasnormalize(las_clip, knnidw(k=10,p=2))
        las_norm_veg=lasfilter(las_norm,Classification==1L)
        
        # Extract horizontal variability and proportion related metrics
        
        dsm = grid_metrics(las_norm,~quantile(Z, 0.90),res=1)
        proj4string(dsm) <- CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs")
        
        dtm=grid_metrics(las_clip,~min(Z[Classification==2L]),res=1)
        proj4string(dtm) <- CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs")
        
        dim = grid_metrics(las_norm,~mean(Intensity),res=1)
        proj4string(dim) <- CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs")
        
        dim_b1 = grid_metrics(las_norm,~mean(Intensity[Z<1]),res=1)
        proj4string(dim_b1) <- CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs")
        
        dsm_b1 = grid_metrics(las_norm,~mean(Z[Z<1]),res=1)
        proj4string(dsm_b1) <- CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs")
        
        ## 
        tri_dsm=terrain(dsm,opt="roughness",neighbors=8)
        tri_dsm_b1=terrain(dsm_b1,opt="roughness",neighbors=8)
        
        var_dsm=focal(dsm,w=matrix(1,3,3), fun=var, pad=TRUE,na.rm = TRUE)
        var_dsm_b1=focal(dsm_b1,w=matrix(1,3,3), fun=var, pad=TRUE,na.rm = TRUE)
        
        var_dim=focal(dim,w=matrix(1,3,3), fun=var, pad=TRUE,na.rm = TRUE)
        var_dim_b1=focal(dim_b1,w=matrix(1,3,3), fun=var, pad=TRUE,na.rm = TRUE)
        
        dtm_slope=terrain(dtm,opt="slope",neighbors=4,unit='degrees')
        dtm_asp=terrain(dtm,opt="aspect",neighbors=4,unit='degrees')
        
        butterflysp_df_sel_sel=butterflysp_df_sel[butterflysp_df_sel$Tr_sec==butterflysp_df_sel$Tr_sec[j],]
        
        coordinates(butterflysp_df_sel_sel)=~x+y
        proj4string(butterflysp_df_sel_sel)<- CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs")
        
        tri_dsm_med <- extract(tri_dsm, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=median)
        tri_dsm_med_b1 <- extract(tri_dsm_b1, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=median)
        
        var_dsm_med <- extract(var_dsm, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=median)
        var_dsm_med_b1 <- extract(var_dsm_b1, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=median)
        
        var_dim_med <- extract(var_dim, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=median)
        var_dim_med_b1 <- extract(var_dim_b1, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=median)
        
        dtm_slope_med <- extract(dtm_slope, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=median)
        dtm_asp_med <- extract(dtm_asp, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=median)
        
        tri_dsm_mean <- extract(tri_dsm, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=mean)
        tri_dsm_mean_b1 <- extract(tri_dsm_b1, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=mean)
        
        var_dsm_mean <- extract(var_dsm, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=mean)
        var_dsm_mean_b1 <- extract(var_dsm_b1, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=mean)
        
        var_dim_mean <- extract(var_dim, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=mean)
        var_dim_mean_b1 <- extract(var_dim_b1, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=mean)
        
        dtm_slope_mean <- extract(dtm_slope, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=mean)
        dtm_asp_mean <- extract(dtm_asp, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=mean)
        
        tri_dsm_sd <- extract(tri_dsm, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=sd)
        tri_dsm_sd_b1 <- extract(tri_dsm_b1, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=sd)
        
        var_dsm_sd <- extract(var_dsm, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=sd)
        var_dsm_sd_b1 <- extract(var_dsm_b1, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=sd)
        
        var_dim_sd <- extract(var_dim, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=sd)
        var_dim_sd_b1 <- extract(var_dim_b1, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=sd)
        
        dtm_slope_sd <- extract(dtm_slope, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=sd)
        dtm_asp_sd <- extract(dtm_asp, butterflysp_df_sel_sel, weights=FALSE,buffer = 100,fun=sd)
        
        height_class=reclassify(dsm, c(-Inf,5,1, 5,Inf,2))
        height_class[height_class == 1] <- NA
        
        propofhighveg <- extract(height_class,butterflysp_df_sel_sel,buffer = 100,fun=length)
        
        if (is.na(propofhighveg)) {propofhighveg=0} else {propofhighveg=propofhighveg}
        
        # Export
        
        newline <- data.frame(t(c(Transect=i,Transect_sec=paste(butterflysp_df_sel$Tr_sec[j],sep=""),
                                  rough_dsm_med=tri_dsm_med,
                                  rough_dsm_med_b1=tri_dsm_med_b1,
                                  var_dsm_med=var_dsm_med,
                                  var_dsm_med_b1=var_dsm_med_b1,
                                  var_dim_med=var_dim_med,
                                  var_dim_med_b1=var_dim_med_b1,
                                  dtm_slope_med=dtm_slope_med,
                                  dtm_asp_med=dtm_asp_med,
                                  rough_dsm_mean=tri_dsm_mean,
                                  rough_dsm_mean_b1=tri_dsm_mean_b1,
                                  var_dsm_mean=var_dsm_mean,
                                  var_dsm_mean_b1=var_dsm_mean_b1,
                                  var_dim_mean=var_dim_mean,
                                  var_dim_mean_b1=var_dim_mean_b1,
                                  dtm_slope_mean=dtm_slope_mean,
                                  dtm_asp_mean=dtm_asp_mean,
                                  rough_dsm_sd=tri_dsm_sd,
                                  rough_dsm_sd_b1=tri_dsm_sd_b1,
                                  var_dsm_sd=var_dsm_sd,
                                  var_dsm_sd_b1=var_dsm_sd_b1,
                                  var_dim_sd=var_dim_sd,
                                  var_dim_med_b1=var_dim_med_b1,
                                  dtm_slope_sd=dtm_slope_sd,
                                  dtm_asp_sd=dtm_asp_sd,
                                  propofhighveg=propofhighveg)))
        
        dpcloudfea_exp_df <- rbind(dpcloudfea_exp_df, newline)
        
      }
    }
  }
}

write.csv(dpcloudfea_exp_df,"Butterfly_lidarmetrics_100m_v3.csv")
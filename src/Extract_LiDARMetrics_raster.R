"
@author: Zsofia Koma, UvA
Aim: LiDAR metrics extraction for the butterfly project
"
library(dplyr)
library(lidR)
library(rgdal)
library(raster)

library(e1071)
library(snow)

shannon = function(Z)
{
  
  p1=(length(Z[(Z>0 & Z<2.5)])/length(Z))*100
  p2=(length(Z[(Z>2.5 & Z<5)])/length(Z))*100
  p3=(length(Z[(Z>5 & Z<7.5)])/length(Z))*100
  p4=(length(Z[(Z>7.5 & Z<10)])/length(Z))*100
  p5=(length(Z[(Z>10 & Z<12.5)])/length(Z))*100
  p6=(length(Z[(Z>12.5 & Z<15)])/length(Z))*100
  p7=(length(Z[(Z>15 & Z<17.5)])/length(Z))*100
  p8=(length(Z[(Z>17.5 & Z<20)])/length(Z))*100
  p9=(length(Z[(Z>20)])/length(Z))*100
  
  v=c(p1,p2,p3,p4,p5,p6,p7,p8,p9)
  
  p <- table(v)
  p <- p/sum(p)
  shannon=sum(-p*log(p))
  
  return(shannon)
}

# Set working dirctory
workingdirectory="D:/Sync/_Amsterdam/08_coauthor_MScProjects/Reinier/lidarmetrics_calc/"
#workingdirectory="D:/Reinier/Reinier_output/"
setwd(workingdirectory)

butterflyspfile="Melitaea_AnalysisData.shp"

#Import
butterflysp = readOGR(dsn=butterflyspfile)
butterflysp_df=butterflysp@data

butterflysp_df_gr <- butterflysp_df %>%
  group_by(Transect) %>%
  summarise(nofobs = length(Transect))

# Direct point cloud based metrics
#Transect <- as.vector(butterflysp_df_gr$Transect)
#Transect= as.numeric(Transect)
Transect=c(123,362)

for (i in Transect) {
  print(i)
  
  las=readLAS(paste("Transect_",i,".laz",sep=""))
  las_norm=lasnormalize(las, knnidw(k=10,p=2))
  writeLAS(las_norm,paste("Transect_",i,"_norm.laz",sep=""))
  
  las_norm_veg=lasfilter(las_norm,Classification==1L)
  
  nofret_pheightlay_b02=grid_metrics(las_norm,~(length(Z[(Classification==1L & Z<0.2)])/length(Z))*100,res=5)
  nofret_pheightlay_02_05=grid_metrics(las_norm,~(length(Z[(Classification==1L & Z>0.2 & Z<0.5)])/length(Z))*100,res=5)
  nofret_pheightlay_05_1=grid_metrics(las_norm,~(length(Z[(Classification==1L & Z>0.5 & Z<1)])/length(Z))*100,res=5)
  nofret_pheightlay_1_2=grid_metrics(las_norm,~(length(Z[(Classification==1L & Z>1 & Z<2)])/length(Z))*100,res=5)
  nofret_pheightlay_2_5=grid_metrics(las_norm,~(length(Z[(Classification==1L & Z>2 & Z<5)])/length(Z))*100,res=5)
  nofret_pheightlay_5_10=grid_metrics(las_norm,~(length(Z[(Classification==1L & Z>5 & Z<10)])/length(Z))*100,res=5)
  nofret_pheightlay_10_20=grid_metrics(las_norm,~(length(Z[(Classification==1L & Z>10 & Z<20)])/length(Z))*100,res=5)
  nofret_pheightlay_a20=grid_metrics(las_norm,~(length(Z[(Classification==1L & Z>20)])/length(Z))*100,res=5)
  
  zmean=grid_metrics(las_norm_veg,~mean(Z),res=5)
  z090quantile = grid_metrics(las_norm_veg,~quantile(Z, 0.90),res=5)
  proj4string(z090quantile) <- CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs")
  
  zmean_undst=grid_metrics(las_norm_veg,~mean(Z[Z<5]),res=5)
  echovar=grid_metrics(las_norm_veg,~var(ReturnNumber),res=5)
  int_mean=grid_metrics(las_norm_veg,~mean(Intensity),res=5)
  int_sd=grid_metrics(las_norm_veg,~sd(Intensity),res=5)
  
  pulsepen = grid_metrics(las_norm,~(length(Z[(Classification==2L)])/length(Z))*100,res=5)
  
  zkurto = grid_metrics(las_norm_veg,~kurtosis(Z),res=5)
  zsd = grid_metrics(las_norm_veg,~sd(Z),res=5)
  z025quantile = grid_metrics(las_norm_veg,~quantile(Z, 0.25),res=5)
  z050quantile = grid_metrics(las_norm_veg,~quantile(Z, 0.50),res=5)
  z075quantile = grid_metrics(las_norm_veg,~quantile(Z, 0.75),res=5)
  
  shannon=grid_metrics(las_norm_veg,~shannon(Z),res=5)
  
  dtm=grid_metrics(las,~min(Z[Classification==2L]),res=5)
  proj4string(dtm) <- CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs")
  
  tri_dsm=terrain(z090quantile,opt="TRI",neighbors=4)
  
  beginCluster(18)
  var_dsm=clusterR(z090quantile, focal, args=list(w=matrix(1,3,3), fun=var, pad=TRUE,na.rm = TRUE))
  endCluster()
  
  names(var_dsm) <- "var_dsm"
  
  dtm_slope=terrain(dtm,opt="slope",neighbors=4,unit='degrees')
  dtm_asp=terrain(dtm,opt="aspect",neighbors=4,unit='degrees')
  
  metrics=stack(nofret_pheightlay_b02,nofret_pheightlay_02_05,nofret_pheightlay_05_1,nofret_pheightlay_1_2,nofret_pheightlay_2_5,nofret_pheightlay_5_10,nofret_pheightlay_10_20,nofret_pheightlay_a20,
                zmean,z090quantile,zmean_undst,echovar,int_mean,int_sd,pulsepen,zkurto,zsd,z025quantile,z050quantile,z075quantile,shannon,dtm,tri_dsm,var_dsm,dtm_slope,dtm_asp)
  
  names(metrics)<- c("nofret_pheightlay_b02","nofret_pheightlay_02_05","nofret_pheightlay_05_1","nofret_pheightlay_1_2","nofret_pheightlay_2_5","nofret_pheightlay_5_10","nofret_pheightlay_10_20","nofret_pheightlay_a20","
                zmean","z090quantile","zmean_undst","echovar","int_mean","int_sd","pulsepen","zkurto","zsd","z025quantile","z050quantile","z075quantile","shannon","dtm","tri_dsm","var_dsm","dtm_slope","dtm_asp")
  
  proj4string(metrics) <- CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs")
  
  writeRaster(metrics,paste("Transect_",i,"_rastermetrics.grd",sep=""),overwrite=TRUE)
  
}
"
@author: Zsofia Koma, UvA
Aim: overlay butterfly transect with LiDAR 
"

# Import libraries
library("rgdal")
library("raster")
library("rgeos")
library("spatialEco")
library("dplyr")

library("ggplot2")
library("gridExtra")

# Set global variables
#full_path="C:/Koma/Sync/_Amsterdam/08_coauthor_MScProjects/Reinier/datapreprocess/"
full_path="D:/Sync/_Amsterdam/08_coauthor_MScProjects/Reinier/datapreprocess/"

transectfile="Transects_0712_Complete.csv"
ahnfile="ahn2.shp"

setwd(full_path)

# Import 

ahn2 = readOGR(dsn=ahnfile)
transect=read.csv(file=transectfile,header=TRUE,sep=";")

# Convert transect to point shp

transect$X=transect$x
transect$Y=transect$y

transect_shp=transect[c("X","Y","x","y","Tr_sec","Transect","Section")]
coordinates(transect_shp)=~X+Y
proj4string(transect_shp)<- CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs")

#rgdal::writeOGR(transect_shp, '.', "transects", 'ESRI Shapefile',overwrite_layer = TRUE)

# Transect points within ahn2

intersect_points_wpoly = spatialEco::point.in.poly(transect_shp, ahn2)
intersect_points_wpoly_df=intersect_points_wpoly@data

write.table(intersect_points_wpoly_df, file = "intersect_points_wpoly_df.csv",row.names=FALSE,col.names=TRUE,sep=",")

per_ahn2 <- intersect_points_wpoly_df %>%
  group_by(bladnr) %>%
  summarise(nofobs = length(bladnr))

## Create one polygon per Transect by extent

transectgroups <- transect %>%
  group_by(Transect) %>%
  summarise(nofobs = length(Transect))

boundary_pertransect=data.frame(Transect=integer(nrow(transectgroups)),xmin=double(nrow(transectgroups)),xmax=double(nrow(transectgroups)),ymin=double(nrow(transectgroups)),ymax=double(nrow(transectgroups)))

for (i in seq_len(nrow(transectgroups))) {
  print(transectgroups$Transect[i])
  
  transect_shp_sel=transect_shp[transect_shp@data$Transect==transectgroups$Transect[i],]
  union=extent(transect_shp_sel)
  
  boundary_pertransect$Transect[i]<-transectgroups$Transect[i]
  boundary_pertransect$xmin[i]<-union@xmin
  boundary_pertransect$xmax[i]<-union@xmax
  boundary_pertransect$ymin[i]<-union@ymin
  boundary_pertransect$ymax[i]<-union@ymax
  
} 

write.table(boundary_pertransect, file = "boundaries_pertransects.csv",row.names=FALSE,col.names=TRUE,sep=",")

# Create a polygon

forpoly=data.frame(Transect=integer(nrow(boundary_pertransect)),wkt_t=character(nrow(boundary_pertransect)),stringsAsFactors=FALSE)

for (i in seq_len(nrow(boundary_pertransect))) {
  print(boundary_pertransect$Transect[i])
  
  boundary_pertransect_sel=boundary_pertransect[boundary_pertransect$Transect==boundary_pertransect$Transect[i],]
  
  boundary_pertransect_sel$xmin=boundary_pertransect_sel$xmin-1000
  boundary_pertransect_sel$xmax=boundary_pertransect_sel$xmax+1000
  boundary_pertransect_sel$ymin=boundary_pertransect_sel$ymin-1000
  boundary_pertransect_sel$ymax=boundary_pertransect_sel$ymax+1000
  
  wkt_astext=paste("POLYGON((",boundary_pertransect_sel$xmin," ",boundary_pertransect_sel$ymin,",",boundary_pertransect_sel$xmin," ",boundary_pertransect_sel$ymax,",", boundary_pertransect_sel$xmax," ",boundary_pertransect_sel$ymax,",",
                   boundary_pertransect_sel$xmax," ", boundary_pertransect_sel$ymin,",",boundary_pertransect_sel$xmin," ",boundary_pertransect_sel$ymin,"))",sep="")
  
  forpoly$Transect[i]<-boundary_pertransect$Transect[i]
  forpoly$wkt_t[i]<-as.character(wkt_astext)
  
}

write.table(forpoly, file = "boundaries_pertransects_wkt_1000.csv",row.names=FALSE,col.names=TRUE,sep=",")


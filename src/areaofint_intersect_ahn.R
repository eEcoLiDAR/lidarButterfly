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
library("maptools")

library("ggplot2")
library("gridExtra")

# Set global variables
#full_path="C:/Koma/Sync/_Amsterdam/08_coauthor_MScProjects/Reinier/datapreprocess/"
full_path="D:/Sync/_Amsterdam/08_coauthor_MScProjects/Reinier/datapreprocess/"

areaofintfile="bound_1km_ptransect.shp"
ahnfile="ahn2.shp"

setwd(full_path)

# Import 

ahn2 = readOGR(dsn=ahnfile)
areaofint=readOGR(dsn=areaofintfile)

# Union of the polygons - do not extract the same km square several times
poly_union <- unionSpatialPolygons(areaofint,rep(1, length(areaofint)))
plot(poly_union)

poly_union <- as(poly_union, "SpatialPolygonsDataFrame")

poly_union_m <- disaggregate(poly_union)
poly_union_m@data$gen_id <- seq(1,length(poly_union_m@data$dummy))

writeOGR(obj=poly_union_m, dsn=".", layer="transect_poly_union", driver="ESRI Shapefile",overwrite_layer = TRUE)

# Intersection
intersected=raster::intersect(ahn2,areaofint)
intersected_df=intersected@data

req_ahn2 <- intersected_df %>%
  group_by(bladnr) %>%
  summarise(nofobs = length(bladnr))

write.table(req_ahn2, file = "req_ahn2_for1km.csv",row.names=FALSE,col.names=TRUE,sep=",")
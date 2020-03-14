"
@author: Zsofia Koma, UvA
Aim: overlay butterfly transect with LiDAR (AHN3)

put into QGIS wkt csv and transfer into shp
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

library(pingr)

# Set global variables
full_path="C:/Koma/Sync/_Amsterdam/08_coauthor_MScProjects/Reinier/dataprocess2/"
#full_path="D:/Sync/_Amsterdam/08_coauthor_MScProjects/Reinier/datapreprocess/"

areaofintfile="C:/Koma/Sync/_Amsterdam/08_coauthor_MScProjects/Reinier/dataprocess2/input/req_transect_groningen_250buffer.shp"
ahnfile="C:/Koma/Sync/_Amsterdam/08_coauthor_MScProjects/Reinier/dataprocess2/input/ahn3.shp"

setwd(full_path)

# Import 

ahn2 = readOGR(dsn=ahnfile)
areaofint=readOGR(dsn=areaofintfile)

# Intersection
intersected=raster::intersect(ahn2,areaofint)
intersected_df=intersected@data

req_ahn2 <- intersected_df %>%
  group_by(bladnr) %>%
  summarise(nofobs = length(bladnr))

req_ahn2$bladnr_up=toupper(req_ahn2$bladnr)

req_ahn2$list=paste("https://geodata.nationaalgeoregister.nl/ahn3/extract/ahn3_laz/C_",req_ahn2$bladnr_up,".LAZ",sep="")

write.table(req_ahn2, file = "req_ahn3_butterfly.csv",row.names=FALSE,col.names=TRUE,sep=",")

#Prepare data for downloading with curl

write.table(req_ahn2$list, file = "ahn3list.txt", append = FALSE, quote = FALSE, sep = "", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double")) 

"
@author: Zsofia Koma, UvA
Aim: Extract area of interest from AHN3 data
"
library("lidR")
library("rgdal")

# Set working dirctory
workingdirectory="D:/Reinier/"
#workingdirectory="C:/Koma/Sync/_Amsterdam/08_coauthor_MScProjects/Reinier/datapreprocess/ahn3/"
setwd(workingdirectory)

#Import 
areaofintfile="req_transect_250buffer.shp"
areaofint=readOGR(dsn=areaofintfile)

# Extract
ctg = catalog(workingdirectory)

for (i in seq(1,length(areaofint$Transect),1)){ 
  print(areaofint$Transect[i]) 
  subset = lasclip(ctg, areaofint[areaofint$Transect==areaofint$Transect[i],])
  
  if (subset@header@PHB[["Number of point records"]]>0) {
    writeLAS(subset,paste("C:/Reinier_output/Transect_",areaofint$Transect[i],".laz",sep=""))
  }
}
"
@author: Zsofia Koma, UvA
Aim: Extract area of interest from AHN3 data

before this you need to download the ahn tiles from the ahn3list using bash
"
#while read p; do curl -o ${p:62:73} "${p%?}";done < ahn3list.txt

library("lidR")
library("rgdal")

# Set working dirctory
workingdirectory="D:/Reinier/Groningen/"
setwd(workingdirectory)

#Import 
areaofintfile="req_transect_groningen_250buffer.shp"
areaofint=readOGR(dsn=areaofintfile)

# Extract
ctg = catalog(workingdirectory)

for (i in seq(1,length(areaofint$Transect),1)){ 
  print(areaofint$Transect[i]) 
  subset = lasclip(ctg, areaofint[areaofint$Transect==areaofint$Transect[i],])
  
  if (subset@header@PHB[["Number of point records"]]>0) {
    writeLAS(subset,paste("D:/Reinier/Groningen/Transect_",areaofint$Transect[i],".laz",sep=""))
  }
}
"
@author: Zsofia Koma, UvA
Aim: Extract area of interest from AHN2 data
"
library("lidR")
library("rgdal")

# Set working dirctory
workingdirectory="D:/Reinier/"
#workingdirectory="D:/Sync/_Amsterdam/08_coauthor_MScProjects/Reinier/datapreprocess/"
setwd(workingdirectory)

#Import 
areaofintfile="transect_poly_union_1kmsquares.shp"
areaofint=readOGR(dsn=areaofintfile)

ctg = catalog(workingdirectory)

for (i in areaofint@data[["id"]]){ 
  print(i) 
  subset = lasclip(ctg, areaofint[areaofint$gen_id==i,])
  
  if (length(subset)>0) {
    for (g in seq(from=1,to=length(subset))) {
      print(g)
      writeLAS(subset[[g]],paste("GenID_",i,"_",g,".laz",sep=""))
    }
  }
}

#i=252
#if (length(subset)>0) {
  #for (g in seq(from=1,to=length(subset))) {
    #print(g)
    #writeLAS(subset[[g]],paste("GenID_",i,"_",g,".laz",sep=""))
    #}
#}

# Clipping
#for (i in seq(from=1,to=max(areaofint$gen_id))){ 
  #print(i)
  
  #subset = lasclip(ctg, areaofint[areaofint$gen_id==i,])
  
  #if (subset@header@PHB[["Number of point records"]]>0) {
    #writeLAS(subset,paste("GenID_",i,".laz",sep=""))
  #}
  
#}
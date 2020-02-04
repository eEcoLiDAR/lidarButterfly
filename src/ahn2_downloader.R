"
@author: Zsofia Koma, UvA
Aim: This script download the required AHN2 data from the pdok website
"

# Set working directory
#full_path="C:/Koma/Sync/_Amsterdam/08_coauthor_MScProjects/Reinier/datapreprocess/"
#full_path="D:/Sync/_Amsterdam/08_coauthor_MScProjects/Reinier/datapreprocess/"
full_path="D:/Reinier/"

reqahnfile="req_ahn2_for1km_v2.csv"

setwd(full_path)

# Import
reqahn=read.csv(file=reqahnfile,header=TRUE,sep=",")

# Set filenames and dwnload and unzip the required dataset

for (tile in reqahn$bladnr){
  print(tile)
  
  download.file(paste("http://geodata.nationaalgeoregister.nl/ahn2/extract/ahn2_","gefilterd/g",tile,".laz.zip",sep=""),
                paste("g",tile,".laz.zip",sep=""))
  download.file(paste("http://geodata.nationaalgeoregister.nl/ahn2/extract/ahn2_","uitgefilterd/u",tile,".laz.zip",sep=""),
                paste("u",tile,".laz.zip",sep=""))
  
  unzip(paste("g",tile,".laz.zip",sep=""))
  unzip(paste("u",tile,".laz.zip",sep=""))
}

zipped <- dir(path=full_path, pattern=".laz.zip")
file.remove(zipped)
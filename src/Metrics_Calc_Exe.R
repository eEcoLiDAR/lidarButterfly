library("lidR")
library("rgdal")

#source("C:/Users/Gebruiker/Documents/R/lidR_butterflies/Metrics_Functions.R")
source("C:/Koma/Github/komazsofi/myPhD_escience_analysis/RV_MSc/Metrics_Functions.R")

# Set Working directory

#workingdirectory="C:/Users/Gebruiker/Documents/R/lidR_butterflies"
workingdirectory="C:/Koma/Paper2/Reinier"
setwd(workingdirectory)

cores=1 #CPU on computer: max-1 (to avoid freezing) #server = 18
chunksize=100
buffer=2.5 # overlap between tiles for neighbourhood calculations e.g. Minimum same as resolution.
resolution=2.5 #resolution of the raster

rasterOptions(maxmemory = 200000000000) #avoid freezing due to reaching max. raster size

# clipping
ctg_clip <- catalog(workingdirectory)

clipped=lasclipCircle(ctg_clip,169844,440463.5,25)
writeLAS(clipped,paste(workingdirectory,"/clipped/Transect_x.laz",sep=""))

ctg <- catalog(paste(workingdirectory,"/clipped/Transect_x.laz",sep=""))
#plot(ctg)

# normalizing
opt_chunk_buffer(ctg) <- buffer
opt_chunk_size(ctg) <- chunksize
opt_cores(ctg) <- cores
opt_output_files(ctg) <- paste(workingdirectory,"/clipped/normalized/{XLEFT}_{YBOTTOM}",sep="")

normalizedctg = lasnormalize(ctg, knnidw(k=20,p=2))

# metrics calculation
opt_output_files(normalizedctg) <- paste(workingdirectory,"/clipped/normalized/{XLEFT}_{YBOTTOM}_height",sep="")
height = grid_metrics(normalizedctg, Heightmetrics(Z),res = resolution)
plot(height)

opt_output_files(normalizedctg) <- paste(workingdirectory,"/clipped/normalized/{XLEFT}_{YBOTTOM}_ampl",sep="")
ampl = grid_metrics(normalizedctg, Amplmetrics(Intensity),res = resolution)
plot(ampl) 
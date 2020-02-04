# Extraction of LidAR points for the lesson

library(lidR)

workdirectory=setwd("C:/Koma/")

catalog=catalog(workdirectory)

clipped=lasclipCircle(catalog,186438,452574,200)
clipped_2=lasclipCircle(catalog,204570,522177,200)
clipped_3=lasclipCircle(catalog,204200,457570,200)
clipped_4=lasclipCircle(catalog,95533,484706,200)

writeLAS(clipped,"362_18.laz")
writeLAS(clipped_2,"460_9.laz")
writeLAS(clipped_3,"375_5.laz")
writeLAS(clipped_4,"469_12.laz")
#### Nitrogen Level fix
Drive<-"C:/Users/zjrobbin/Documents/GitHub/Project-Southern-Appalachians-2018/Models/LANDIS_Sapps_Active_v1_2/"
MRSOM1surfN<-raster(paste0(Drive,"MRSOM1surfNmap.tif"))
MRSOM1soilN<-raster(paste0(Drive,"MRSOM1surfNmap.tif"))
plot(MRSOM1surfN,zlim=c(0,10))

##Take MRSOM1surfN to 2/3
MRSOM1surfN_updated<-MRSOM1surfN*.1
plot(MRSOM1surfN_updated,zlim=c(0,10))
## Take MRSom1Soil to Half. 
plot(MRSOM1soilN, zlim=c(0,10))
MRSOM1soilN_updated<-MRSOM1soilN*.1
plot(MRSOM1soilN_updated,zlim=c(0,10))

writeRaster(MRSOM1soilN_updated,paste0(Drive,"MRSOM1soilN_updated.tif"),overwrite=T)
writeRaster(MRSOM1surfN_updated,paste0(Drive,"MRSOM1surfN_updated.tif"),overwrite=T)

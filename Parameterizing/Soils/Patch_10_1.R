library(raster)
Drive<-'C:/Users/zacha/Documents/GitHub/Project-Southern-Appalachians-2018/Parameterizing/Soils/'
Surf<-raster::raster(paste0(Drive,"MRSOM1surfCmap.tif"))
plot(Surf,zlim=c(0,100))
SurfN<-Surf/10
plot(SurfN,zlim=c(0,10))

Soil<-raster::raster(paste0(Drive,"MRSOM1SoilCmap.tif"))
SoilN<-Soil/10
plot(SoilN,zlim=c(0,12))

Soil2<-raster::raster(paste0(Drive,"MRSOM2Cmap.tif"))
Soil2N<-Soil2/17.5




writeRaster(SoilN,paste0(Drive,"SOM1soilNmap_v1_1.TIF"),overwrite=T)
writeRaster(SurfN,paste0(Drive,'SOM1surfNmap_v1_1.TIF'),overwrite=T)
writeRaster(Soil2N,paste0(Drive,'SOM2Nmap_v1_1.TIF'),overwrite=T)


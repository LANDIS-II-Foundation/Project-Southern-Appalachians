library(raster)
Drive<-'C:/Users/zacha/Documents/GitHub/Project-Southern-Appalachians-2018/Parameterizing/Soils/'
Surf<-raster::raster(paste0(Drive,"MRSOM1surfCmap.tif"))
Surf2C<-Surf*.07                     
plot(Surf2,zlim=c(0,5))
SurfN<-Surf2/100
plot(SurfN,zlim=c(0,.1))

Soil<-raster::raster(paste0(Drive,"MRSOM1SoilCmap.tif"))
SoilN<-Soil/8
plot(SoilN,zlim=c(0,12))

Soil2<-raster::raster(paste0(Drive,"MRSOM2SoilCmap.tif"))





writeRaster(SoilN,paste0(Drive,"SOM1soilNmap_v1_1.tif"))
writeRaster(SurfN,paste0(Drive,'SOM1surfNmap_v1_1.tif'))
writeRaster(Surf2C,paste0(Drive,'SOM1surfCmap_v1_1.tif'))

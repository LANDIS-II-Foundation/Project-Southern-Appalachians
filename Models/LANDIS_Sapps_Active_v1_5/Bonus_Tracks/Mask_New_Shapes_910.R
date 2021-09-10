library(sf)
library(raster)
ecos<-raster("C:/Users/zacha/Documents/GitHub/Project-Southern-Appalachians-2018/Models/LANDIS_Sapps_Active_v1_5/Ecos11_NLCD.tif")
plot(ecos)
Mask<-sf::st_read("C:/Users/zacha/Desktop/Sapps_Mortality_Clean/Georgia.shp")
Mask2<-st_transform(Mask,crs=projection(ecos))
cut<-mask(ecos,Mask2) 
plot(cut)
writeRaster(cut,"C:/Users/zacha/Documents/GitHub/Project-Southern-Appalachians-2018/Models/LANDIS_Sapps_Active_v1_3/Georgia_Ecos.tif",
            datatype='INT4S',overwrite=T)


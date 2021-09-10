
ecos<-raster("C:/Users/zacha/Documents/GitHub/Project-Southern-Appalachians-2018/Models/LANDIS_Sapps_Active_v1_3/11_Ecoregions.tif")
plot(ecos)
library(sf)
#mask<-
Mask<-sf::st_read("C:/Users/zacha/Desktop/Sapps_Mortality_Clean/Georgia.shp")

Mask2<-st_transform(Mask,crs=projection(ecos))
cut<-mask(ecos,Mask2) 
#dtypes(ecos)
#dataType(ecos)
cut[is.na(cut$X11_Ecoregions)]<-1
projection(cut)<-projection(ecos)
length(cut$X11_Ecoregions[cut$X11_Ecoregions!=1])
99688/1400700
1400700*(6.25)

#unique(cut$X11_Ecoregions)
plot(cut)
#plot(c
writeRaster(cut,"C:/Users/zacha/Documents/GitHub/Project-Southern-Appalachians-2018/Models/LANDIS_Sapps_Active_v1_3/Georgia_Ecos.tif",
            datatype='INT4S',overwrite=T)

?writeRaster

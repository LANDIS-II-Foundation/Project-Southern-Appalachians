### This script rasterizes the human accidental ignitions from Karen Short's database. Within 1k cells, it counts all of the historic ignitions. 
### Cells with 0 count yield NA, but are forced to zero. The 1km raster is resampled to reflect ~260m cells of the study area, then scaled 0 to 1. 
### Since values of 0 are not desirable for ignition probabilities (because ignitions would never be ignited in the model), we bump the zeros up to 0.05
### because 1/20 years (our period of record), yields that many significant figures. 

#packages for efficiently handling large tabular and spatial data
library(raster)
library(data.table)
library(rgdal)
library(sf)

setwd("C:\\Users\\thebrain\\Dropbox\\DEAL_lab\\SCRPPLE")


#read in short ignition data
fire<-st_read("C:\\Users\\thebrain\\Dropbox\\Firechasers","fire_nalcc")
fire_dt<-as.data.table(fire)
fires<-fire_dt

#all fires types listed, need to subset down to causes that reflect "human accidental"
fire_types<-as.character(unique(fires$STAT_CAU_1)) #all causes
non_humanacc <- c("Lightning","Miscellaneous",'Missing/Undefined') #make list of those to remove
human_acc<-setdiff(fire_types, non_humanacc) #difference between the two lists

#in tabular form, make it smaller
subset_fires<-fires[fires$STATE=="GA" | fires$STATE=="TN" | fires$STATE=="SC" | fires$STATE=="NC",]

#query the fires for causes that match the human accidental list 
human_acc<-subset_fires[subset_fires$STAT_CAUSE %in% human_acc,]

#assign lat/longs for making spatial object coords
lat<-human_acc$LATITUDE
long<-human_acc$LONGITUDE

#bind them up
fire_coords<-cbind(long,lat)

#make sure wgs 1984, matching shorts lat/longs
fire_pnts_wgs<-SpatialPoints(fire_coords,proj4string=CRS(as.character("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")))

#clip 4 states worth of ignitions down to the study area extent
aoi<-raster("C:\\Users\\thebrain\\Dropbox\\DEAL_lab\\S_Apps_Project\\Project-Southern-Appalachians-2018\\Models\\LANDIS_Sapps_Active_v1_3\\11_Ecoregions.tif")
#get the specs for the study area maps
utm<-crs(aoi)
crs(aoi)
e <- extent(aoi)

# make a SpatialPolygons object from the study area extent information
p <- as(e, 'SpatialPolygons') 
crs(p)<-proj4string(aoi)


#now transform the points from wgs to utm
fire_utm<-spTransform(fire_pnts_wgs,CRS(proj4string(aoi)))

#now, since both are in correct, metered project - clip!
fires_inextent<-crop(fire_utm,p)

#write out the fires/ignitions as points within our study area so we have it for later
shapefile(fires_inextent,"fire_points_aoi.shp")

###now for rasterizing

#rasterize those points to prespecified raster template
rast<-raster()
crs(rast)<-crs(aoi)
extent(rast)<-extent(aoi)
#res(rast)<-res(aoi)
res(rast)<-1000 #choose to rasterize at 1km  scale because you get higher ignitions counts than at 250m scale 

#rasterize the ignition count to the raster we created just above
ha_surf<-rasterize(fires_inextent,rast,fun="count")

#get rid of NAs
ha_surf[is.na(ha_surf)]<-0

#disaggregate the data without interpolation
agg1000_thendis250<-disaggregate(ha_surf, fact=4)

#scale the data from 0 to 1
ha_scale<-rescale0to1(agg1000_thendis250)

#can't have 0 values, so want to just barely assign the mimnimum value
ha_scale[ha_scale==0]<-0.01

#resample to get to correct study area raster specs
test<-resample(ha_out,aoi,method="bilinear")

#round it and multiple by 100
ha_out<-round(test,3)*100 #multiple by 100 because it's a relative weight

#write it out
writeRaster(ha_out,"HA_ign.tif",overwrite=T)





### everything below used to test best scaling/normlization/resolution


#create function for min/max normalization
# normalize <- function(x) {
#   
#   min <- raster::minValue(x)
#   max <- raster::maxValue(x)
#   
#   return ((x - min) / (max - min))
# }

#use function for min max normalization
# ha_norm<-normalize(ha_surf)
# #now scale those values 0 - 1, which is what we need for the probability raster
# ha_scale<-rescale0to1(ha_surf)



# ha_norm_250<-ha_norm
# ha_scale_250<-ha_scale
# ha_norm_1000<-raster::aggregate(ha_norm,fact=4)
# ha_scale_1000<-raster::aggregate(ha_scale,fact=4)
# ha_norm_5000<-raster::aggregate(ha_norm,fact=20)
# ha_scale_5000<-raster::aggregate(ha_scale,fact=20)


#used these to test different scaling/normalization/resampling

# writeRaster(ha_norm,"ha_norm_500.tif",overwrite=T)
# writeRaster(ha_scale,"ha_scale_500.tif",overwrite=T)
# 
# writeRaster(ha_norm_250,"ha_norm_250.tif",overwrite=T)
# writeRaster(ha_scale_250,"ha_scale_250.tif",overwrite=T)
# 
# 
# 
# writeRaster(ha_norm_1000,"ha_norm_1000.tif",overwrite=T)
# writeRaster(ha_scale_1000,"ha_scale_1000.tif",overwrite=T)
# 
# writeRaster(ha_norm_5000,"ha_norm_5000.tif",overwrite=T)
# writeRaster(ha_scale_5000,"ha_scale_5000.tif",overwrite=T)




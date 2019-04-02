library(raster)
library(rgdal)

#set directory
directory <- setwd("C:/Users/ctgerstl/Documents/LANDIS_Sapps_PracticeLandscape_1_23")

#inport maps
EcoregionMap <- raster(paste(directory, "/MR_DummyEco.tif", sep = ""))
SoilDepthMap <- raster(paste(directory, "/MRSoil_Depth_done.tif", sep = ""))

#Check if projection and Extent are equal
projection(SoilDepthMap) == projection(EcoregionMap)
extent(SoilDepthMap) == extent(EcoregionMap)

#Convert Rasters to Vectors and c-bind them into a dataframe
EcoregionMap.vector <- as.vector(EcoregionMap)
SoilDepthMap.vector <- as.vector(SoilDepthMap)

DataFrame <- as.data.frame(cbind(EcoregionMap.vector, SoilDepthMap.vector))

#Create a new Ecoregion Map vector where the values are set to 2 if Soil depth = 0
DataFrame$NewEcoregionMap <- ifelse(DataFrame$SoilDepthMap.vector<=0, yes=2, no=DataFrame$EcoregionMap.vector)
DataFrame$NewEcoregionMap <- ifelse(is.na(DataFrame$NewEcoregionMap), yes=2, no=DataFrame$NewEcoregionMap)

#Isolate the new values into a vector
NewEcoregionMap.vector <- as.integer(DataFrame$NewEcoregionMap)

#Create an empty raster with the same dimensions, resolution, and crs as the other rasters
Extent <- extent(EcoregionMap)
Resolution <- res(EcoregionMap)
Projection <- projection(EcoregionMap)

NewEcoregionMap.raster <- raster(crs=Projection, ext=Extent, res=Resolution)

#Add in values to the new raster layer and check to make sure it's correct by plotting it
NewEcoregionMap.raster[] <- NewEcoregionMap.vector
plot(NewEcoregionMap.raster)

#Export new raster as a .tif file for use in LANDIS-II runs
writeRaster(NewEcoregionMap.raster, filename=(paste(directory, "/MR_DummyEco_ChrisFixed.tif", sep="")), format="GTiff", overwrite=TRUE, datatype="INT2S")

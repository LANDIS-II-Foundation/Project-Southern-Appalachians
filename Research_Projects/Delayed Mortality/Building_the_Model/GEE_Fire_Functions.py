import ee
import pandas
import geopandas as gp 
import folium
import geemap
import pandas as pd
import math 
print(ee.__version__)
ee.Initialize()
###Process DRdNBR
def Cal_DRdNBR(MTBS_Row):
	shapefile =MTBS_Row
	features = []
	for i in range(shapefile.shape[0]):
		geom = shapefile.iloc[i:i+1,:] 
		jsonDict = eval(geom.to_json()) 
		geojsonDict = jsonDict['features'][0] 
		features.append(ee.Feature(geojsonDict)) 

	AOI = ee.FeatureCollection(features)
	#bb_one=MTBS_Row.total_bounds
	#AOI= ee.Geometry.Polygon(bb_one.tolist())
	#Date=(pd.to_datetime())
	#print(Date)

	Month=str(int(MTBS_Row.StartMonth))
	Year=int(MTBS_Row.Year)    
    
		### Setting up Date Framework
	### If a fire happens after Jun we will take the that year (May-Fire) and the year previous as the baseline, if it happens before leaf 
	### out we will take the previous two years (Growing season only May 1 to Oct )
	if(int(Month)>=10):
		YearZeroEnd=pd.to_datetime(str(int(MTBS_Row.Year))+'-'+str(10)+'-'+str(15))
		YearZeroStart=pd.to_datetime(str(int(MTBS_Row.Year))+'-'+str(5)+'-'+str(1))
		YearSubOneEnd=pd.to_datetime(str(Year-1)+'-'+str(10)+'-'+str(15))
		YearSubOneStart=pd.to_datetime(str(Year-1)+'-'+str(5)+'-'+str(1))
		YearOneEnd=pd.to_datetime(str(Year+1)+'-'+str(10)+'-'+str(15))
		YearOneStart=pd.to_datetime(str(Year+1)+'-'+str(5)+'-'+str(15))
		YearTwoEnd=pd.to_datetime(str(Year+2)+'-'+str(10)+'-'+str(15))
		YearTwoStart=pd.to_datetime(str(Year+2)+'-'+str(5)+'-'+str(15))
		YearThreeEnd=pd.to_datetime(str(Year+3)+'-'+str(10)+'-'+str(15))
		YearThreeStart=pd.to_datetime(str(Year+3)+'-'+str(5)+'-'+str(15))
	if(int(Month)>6 and int(Month)<10):
		YearZeroEnd=pd.to_datetime(str(int(MTBS_Row.Year))+'-'+str(int(MTBS_Row.StartMonth))+'-'+str(int(MTBS_Row.StartDay)))
		YearZeroStart=pd.to_datetime(str(int(MTBS_Row.Year))+'-'+str(5)+'-'+str(1))
		YearSubOneEnd=pd.to_datetime(str(Year-1)+'-'+str(10)+'-'+str(15))
		YearSubOneStart=pd.to_datetime(str(Year-1)+'-'+str(5)+'-'+str(1))
		YearOneEnd=pd.to_datetime(str(Year)+'-'+str(10)+'-'+str(15))
		YearOneStart=pd.to_datetime(str(Year)+'-'+str(int(MTBS_Row.StartMonth))+'-'+str(int(MTBS_Row.StartDay)))
		YearTwoEnd=pd.to_datetime(str(Year+1)+'-'+str(10)+'-'+str(15))
		YearTwoStart=pd.to_datetime(str(Year+1)+'-'+str(5)+'-'+str(15))
		YearThreeEnd=pd.to_datetime(str(Year+2)+'-'+str(10)+'-'+str(15))
		YearThreeStart=pd.to_datetime(str(Year+2)+'-'+str(5)+'-'+str(15))    
	if(int(Month)<=6):
		YearSubOneEnd=pd.to_datetime(str(Year-2)+'-'+str(10)+'-'+str(15))
		YearSubOneStart=pd.to_datetime(str(Year-2)+'-'+str(5)+'-'+str(1))
		YearZeroEnd=pd.to_datetime(str(Year-1)+'-'+str(10)+'-'+str(15))
		YearZeroStart=pd.to_datetime(str(Year-1)+'-'+str(5)+'-'+str(1))
		YearOneEnd=pd.to_datetime(str(Year)+'-'+str(10)+'-'+str(15))
		YearOneStart=pd.to_datetime(str(Year)+'-'+str(int(MTBS_Row.StartMonth))+'-'+str(int(MTBS_Row.StartDay)))
		YearTwoEnd=pd.to_datetime(str(Year+2)+'-'+str(10)+'-'+str(15))
		YearTwoStart=pd.to_datetime(str(Year+2)+'-'+str(5)+'-'+str(15))
		YearThreeEnd=pd.to_datetime(str(Year+3)+'-'+str(10)+'-'+str(15))
		YearThreeStart=pd.to_datetime(str(Year+3)+'-'+str(5)+'-'+str(15))
	def maskL8sr(image):
	 # // Bits 3 and 5 are cloud shadow and cloud, respectively.
	  cloudShadowBitMask = 1 << 3;
	  cloudsBitMask = 1 << 5;
	  snowBitMask = 1 << 4;
	  #// Get the pixel QA band.
	  qa = image.select('pixel_qa');
	  #// All flags should be set to zero, indicating clear conditions.
	  mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)and(qa.bitwiseAnd(cloudsBitMask).eq(0))and(qa.bitwiseAnd(snowBitMask).eq(0));
	  ##// Return the masked image, scaled to TOA reflectance, without the QA bands.
	  return image.updateMask(mask).select("B[0-9]*").copyProperties(image, ["system:time_start"]);

	def getQABits(image, start, end, mascara):
		# Compute the bits we need to extract.
		pattern = 0
		for i in range(start,end+1):
			pattern += 2**i
		# Return a single band image of the extracted QA bits, giving the     band a new name.
		return image.select([0], [mascara]).bitwiseAnd(pattern).rightShift(start)
	#A function to mask out cloudy pixels.
	def maskQuality(image):
		# Select the QA band.
		QA = image.select('pixel_qa')
		# Get the internal_cloud_algorithm_flag bit.
		sombra = getQABits(QA,3,3,'cloud_shadow')
		nubes = getQABits(QA,5,5,'cloud')
		#  var cloud_confidence = getQABits(QA,6,7,  'cloud_confidence')
		cirrus_detected = getQABits(QA,9,9,'cirrus_detected')
		#var cirrus_detected2 = getQABits(QA,8,8,  'cirrus_detected2')
		#Return an image masking out cloudy areas.
		return image.updateMask(sombra.eq(0)).updateMask(nubes.eq(0).updateMask(cirrus_detected.eq(0)))

	if Year >=2015:
		Ls8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR");
		
		#def ProcessOut(Ls8)
		### Pull first Prefire Image

		prefireImCol1 = ee.ImageCollection(Ls8
			#// Filter by dates.
			.filterDate(YearZeroStart, YearZeroEnd)
			#// Filter by location.
			.filterBounds(AOI));

		### Pull first Prefire Image
		prefireImCol2 = ee.ImageCollection(Ls8
			#// Filter by dates.
			.filterDate(YearSubOneStart, YearSubOneEnd)
			#// Filter by location.
			.filterBounds(AOI));
		PreFireImages = prefireImCol1.merge(prefireImCol2);

		prefire_CM_ImCol=PreFireImages.map(maskQuality)
		#prefire_CM_ImCol = PreFireImages.map(maskL8sr);
		pre_cm_mos =prefire_CM_ImCol.median().clip(AOI);
		#print(AOI.centroid())

		#### Post Fire 
		postfireImCol1 = ee.ImageCollection(Ls8
			#// Filter by dates.
			.filterDate(YearOneStart, YearOneEnd)
			#// Filter by location.
			.filterBounds(AOI));
		postfireImCol2 = ee.ImageCollection(Ls8
			#// Filter by dates.
			.filterDate(YearTwoStart, YearTwoEnd)
			#// Filter by location.
			.filterBounds(AOI));
		#count = postfireImCol2.size();
		#print('Count2: ', count);
		postfireImCol3 = ee.ImageCollection(Ls8
			#// Filter by dates.
			.filterDate(YearThreeStart, YearThreeEnd)
			#// Filter by location.
			.filterBounds(AOI));
		#count = postfireImCol3.size()
		#print('Count3: ', count);

		postfireImCol1=postfireImCol1.map(maskQuality).median().clip(AOI)
		postfireImCol2=postfireImCol2.map(maskQuality).median().clip(AOI)
		postfireImCol3=postfireImCol3.map(maskQuality).median().clip(AOI)

		PostfireImages=ee.ImageCollection([postfireImCol1, postfireImCol2,postfireImCol3]);

		### This mosaic needs to be fixed (it currently is best quality last to first), ## Needs to be median. 
		def mapNBR(image):
			NBR=image.normalizedDifference(['B5', 'B7'])
			return(NBR)
		### This mosaic needs to be fixed (it currently is best quality last to first), ## Needs to be median. 
		postNBRcol =PostfireImages.map(mapNBR)
		postNBR=postNBRcol.min().clip(AOI)
		preNBR = pre_cm_mos.normalizedDifference(['B5', 'B7']);
		#postNBR = post_cm_mos.normalizedDifference(['B5', 'B7']);
		dNBR_unscaled = preNBR.subtract(postNBR)
		RdNBR=dNBR_unscaled.divide(preNBR.sqrt()).multiply(1000)
	else: 	
		Ls7 = ee.ImageCollection("LANDSAT/LE07/C01/T1_SR")
		#NBR is used to identify burned areas and provide a measure of burn severity. It is calculated as a ratio between the NIR and SWIR values in traditional fashion.
		#(NIR - SWIR) / (NIR + SWIR)
		#In Landsat 4-7, NBR = (Band 4 – Band 7) / (Band 4 + Band 7).

		prefireImCol1 = ee.ImageCollection(Ls7
			#// Filter by dates.
			.filterDate(YearZeroStart, YearZeroEnd)
			#// Filter by location.
			.filterBounds(AOI));


		### Pull first Prefire Image
		prefireImCol2 = ee.ImageCollection(Ls7
			#// Filter by dates.
			.filterDate(YearSubOneStart, YearSubOneEnd)
			#// Filter by location.
			.filterBounds(AOI));
		PreFireImages = prefireImCol1.merge(prefireImCol2);


		prefire_CM_ImCol=PreFireImages.map(maskQuality)
		#prefire_CM_ImCol = PreFireImages.map(maskL8sr);
		pre_cm_mos =prefire_CM_ImCol.median().clip(AOI);
		#print(AOI.centroid())

		#### Post Fire 
		postfireImCol1 = ee.ImageCollection(Ls7
			#// Filter by dates.
			.filterDate(YearOneStart, YearOneEnd)
			#// Filter by location.
			.filterBounds(AOI));
		postfireImCol2 = ee.ImageCollection(Ls7
			#// Filter by dates.
			.filterDate(YearTwoStart, YearTwoEnd)
			#// Filter by location.
			.filterBounds(AOI));
		#count = postfireImCol2.size();
		#print('Count2: ', count);
		postfireImCol3 = ee.ImageCollection(Ls7
			#// Filter by dates.
			.filterDate(YearThreeStart, YearThreeEnd)
			#// Filter by location.
			.filterBounds(AOI));
		#count = postfireImCol3.size()
		#print('Count3: ', count);

		postfireImCol1=postfireImCol1.map(maskQuality).median().clip(AOI)
		postfireImCol2=postfireImCol2.map(maskQuality).median().clip(AOI)
		postfireImCol3=postfireImCol3.map(maskQuality).median().clip(AOI)

		PostfireImages=ee.ImageCollection([postfireImCol1, postfireImCol2,postfireImCol3]);

		
		def mapNBR(image):
			NBR=image.normalizedDifference(['B5', 'B7'])
			return(NBR)
		### This mosaic needs to be fixed (it currently is best quality last to first), ## Needs to be median. 
		postNBRcol =PostfireImages.map(mapNBR)
		postNBR=postNBRcol.min().clip(AOI)
		preNBR = pre_cm_mos.normalizedDifference(['B5', 'B7']);
		#postNBR = post_cm_mos.normalizedDifference(['B5', 'B7']);
		dNBR_unscaled = preNBR.subtract(postNBR)
		RdNBR=dNBR_unscaled.divide(preNBR.sqrt()).multiply(1000)
	return(preNBR,postNBR,RdNBR)
################ Calculating one year RDNBR
def Cal_1yrDRdNBR(MTBS_Row):
	shapefile =MTBS_Row
	features = []
	for i in range(shapefile.shape[0]):
		geom = shapefile.iloc[i:i+1,:] 
		jsonDict = eval(geom.to_json()) 
		geojsonDict = jsonDict['features'][0] 
		features.append(ee.Feature(geojsonDict)) 

	AOI = ee.FeatureCollection(features)
	#bb_one=MTBS_Row.total_bounds
	#AOI= ee.Geometry.Polygon(bb_one.tolist())
	#Date=(pd.to_datetime())
	#print(Date)

	Month=str(int(MTBS_Row.StartMonth))
	Year=int(MTBS_Row.Year)    
    
		### Setting up Date Framework
	### If a fire happens after Jun we will take the that year (May-Fire) and the year previous as the baseline, if it happens before leaf 
	### out we will take the previous two years (Growing season only May 1 to Oct )
	if(int(Month)>=10):
		YearZeroEnd=pd.to_datetime(str(int(MTBS_Row.Year))+'-'+str(10)+'-'+str(15))
		YearZeroStart=pd.to_datetime(str(int(MTBS_Row.Year))+'-'+str(5)+'-'+str(1))
		YearSubOneEnd=pd.to_datetime(str(Year-1)+'-'+str(10)+'-'+str(15))
		YearSubOneStart=pd.to_datetime(str(Year-1)+'-'+str(5)+'-'+str(1))
		YearOneEnd=pd.to_datetime(str(Year+1)+'-'+str(10)+'-'+str(15))
		YearOneStart=pd.to_datetime(str(Year+1)+'-'+str(5)+'-'+str(15))
		
	if(int(Month)>6 and int(Month)<10):
		YearZeroEnd=pd.to_datetime(str(int(MTBS_Row.Year))+'-'+str(int(MTBS_Row.StartMonth))+'-'+str(int(MTBS_Row.StartDay)))
		YearZeroStart=pd.to_datetime(str(int(MTBS_Row.Year))+'-'+str(5)+'-'+str(1))
		YearSubOneEnd=pd.to_datetime(str(Year-1)+'-'+str(10)+'-'+str(15))
		YearSubOneStart=pd.to_datetime(str(Year-1)+'-'+str(5)+'-'+str(1))
		YearOneEnd=pd.to_datetime(str(Year)+'-'+str(10)+'-'+str(15))
		YearOneStart=pd.to_datetime(str(Year)+'-'+str(int(MTBS_Row.StartMonth))+'-'+str(int(MTBS_Row.StartDay)))
		  
	if(int(Month)<=6):
		YearSubOneEnd=pd.to_datetime(str(Year-2)+'-'+str(10)+'-'+str(15))
		YearSubOneStart=pd.to_datetime(str(Year-2)+'-'+str(5)+'-'+str(1))
		YearZeroEnd=pd.to_datetime(str(Year-1)+'-'+str(10)+'-'+str(15))
		YearZeroStart=pd.to_datetime(str(Year-1)+'-'+str(5)+'-'+str(1))
		YearOneEnd=pd.to_datetime(str(Year)+'-'+str(10)+'-'+str(15))
		YearOneStart=pd.to_datetime(str(Year)+'-'+str(int(MTBS_Row.StartMonth))+'-'+str(int(MTBS_Row.StartDay)))
		
	def maskL8sr(image):
	 # // Bits 3 and 5 are cloud shadow and cloud, respectively.
	  cloudShadowBitMask = 1 << 3;
	  cloudsBitMask = 1 << 5;
	  snowBitMask = 1 << 4;
	  #// Get the pixel QA band.
	  qa = image.select('pixel_qa');
	  #// All flags should be set to zero, indicating clear conditions.
	  mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)and(qa.bitwiseAnd(cloudsBitMask).eq(0))and(qa.bitwiseAnd(snowBitMask).eq(0));
	  ##// Return the masked image, scaled to TOA reflectance, without the QA bands.
	  return image.updateMask(mask).select("B[0-9]*").copyProperties(image, ["system:time_start"]);

	def getQABits(image, start, end, mascara):
		# Compute the bits we need to extract.
		pattern = 0
		for i in range(start,end+1):
			pattern += 2**i
		# Return a single band image of the extracted QA bits, giving the     band a new name.
		return image.select([0], [mascara]).bitwiseAnd(pattern).rightShift(start)
	#A function to mask out cloudy pixels.
	def maskQuality(image):
		# Select the QA band.
		QA = image.select('pixel_qa')
		# Get the internal_cloud_algorithm_flag bit.
		sombra = getQABits(QA,3,3,'cloud_shadow')
		nubes = getQABits(QA,5,5,'cloud')
		#  var cloud_confidence = getQABits(QA,6,7,  'cloud_confidence')
		cirrus_detected = getQABits(QA,9,9,'cirrus_detected')
		#var cirrus_detected2 = getQABits(QA,8,8,  'cirrus_detected2')
		#Return an image masking out cloudy areas.
		return image.updateMask(sombra.eq(0)).updateMask(nubes.eq(0).updateMask(cirrus_detected.eq(0)))

	if Year >=2015:
		Ls8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR");
		
		#def ProcessOut(Ls8)
		### Pull first Prefire Image

		prefireImCol1 = ee.ImageCollection(Ls8
			#// Filter by dates.
			.filterDate(YearZeroStart, YearZeroEnd)
			#// Filter by location.
			.filterBounds(AOI));

		### Pull first Prefire Image
		prefireImCol2 = ee.ImageCollection(Ls8
			#// Filter by dates.
			.filterDate(YearSubOneStart, YearSubOneEnd)
			#// Filter by location.
			.filterBounds(AOI));
		PreFireImages = prefireImCol1.merge(prefireImCol2);

		prefire_CM_ImCol=PreFireImages.map(maskQuality)
		#prefire_CM_ImCol = PreFireImages.map(maskL8sr);
		pre_cm_mos =prefire_CM_ImCol.median().clip(AOI);
		#print(AOI.centroid())

		#### Post Fire 
		postfireImCol1 = ee.ImageCollection(Ls8
			#// Filter by dates.
			.filterDate(YearOneStart, YearOneEnd)
			#// Filter by location.
			.filterBounds(AOI));

		#count = postfireImCol3.size()
		#print('Count3: ', count);

		postfireImCol1=postfireImCol1.map(maskQuality).median().clip(AOI)
	

		PostfireImages=ee.ImageCollection([postfireImCol1]);

		def mapNBR(image):
			NBR=image.normalizedDifference(['B5', 'B7'])
			return(NBR)
		### This mosaic needs to be fixed (it currently is best quality last to first), ## Needs to be median. 
		postNBRcol =PostfireImages.map(mapNBR)
		postNBR=postNBRcol.min().clip(AOI)
		preNBR = pre_cm_mos.normalizedDifference(['B5', 'B7']);
		#postNBR = post_cm_mos.normalizedDifference(['B5', 'B7']);
		dNBR_unscaled = preNBR.subtract(postNBR)
		RdNBR=dNBR_unscaled.divide(preNBR.sqrt()).multiply(1000)
	else: 	
		Ls7 = ee.ImageCollection("LANDSAT/LE07/C01/T1_SR")
		#NBR is used to identify burned areas and provide a measure of burn severity. It is calculated as a ratio between the NIR and SWIR values in traditional fashion.
		#(NIR - SWIR) / (NIR + SWIR)
		#In Landsat 4-7, NBR = (Band 4 – Band 7) / (Band 4 + Band 7).

		prefireImCol1 = ee.ImageCollection(Ls7
			#// Filter by dates.
			.filterDate(YearZeroStart, YearZeroEnd)
			#// Filter by location.
			.filterBounds(AOI));


		### Pull first Prefire Image
		prefireImCol2 = ee.ImageCollection(Ls7
			#// Filter by dates.
			.filterDate(YearSubOneStart, YearSubOneEnd)
			#// Filter by location.
			.filterBounds(AOI));
		PreFireImages = prefireImCol1.merge(prefireImCol2);


		prefire_CM_ImCol=PreFireImages.map(maskQuality)
		#prefire_CM_ImCol = PreFireImages.map(maskL8sr);
		pre_cm_mos =prefire_CM_ImCol.median().clip(AOI);
		#print(AOI.centroid())

		#### Post Fire 
		postfireImCol1 = ee.ImageCollection(Ls7
			#// Filter by dates.
			.filterDate(YearOneStart, YearOneEnd)
			#// Filter by location.
			.filterBounds(AOI));
		

		postfireImCol1=postfireImCol1.map(maskQuality).median().clip(AOI)
		

		PostfireImages=ee.ImageCollection([postfireImCol1]);
		
		def mapNBR(image):
			NBR=image.normalizedDifference(['B5', 'B7'])
			return(NBR)
		### This mosaic needs to be fixed (it currently is best quality last to first), ## Needs to be median. 
		postNBRcol =PostfireImages.map(mapNBR)
		postNBR=postNBRcol.min().clip(AOI)
		preNBR = pre_cm_mos.normalizedDifference(['B5', 'B7']);
		#postNBR = post_cm_mos.normalizedDifference(['B5', 'B7']);
		dNBR_unscaled = preNBR.subtract(postNBR)
		RdNBR=dNBR_unscaled.divide(preNBR.sqrt()).multiply(1000)
	return(preNBR,postNBR,RdNBR)



	
	
	
	
##################### variable Processes
def MODIS16(MTBS_Row):
	shapefile =MTBS_Row
	features = []
	for i in range(shapefile.shape[0]):
		geom = shapefile.iloc[i:i+1,:] 
		jsonDict = eval(geom.to_json()) 
		geojsonDict = jsonDict['features'][0] 
		features.append(ee.Feature(geojsonDict)) 

	AOI = ee.FeatureCollection(features)
	#bb_one=MTBS_Row.total_bounds
	#AOI= ee.Geometry.Polygon(bb_one.tolist())
	#Date=(pd.to_datetime())
	#print(Date)

	Month=str(int(MTBS_Row.StartMonth))
	Year=int(MTBS_Row.Year)    
    
	PullDate=pd.to_datetime(str(int(MTBS_Row.Year))+'-'+str(int(MTBS_Row.StartMonth))+'-'+str(int(MTBS_Row.StartDay)))
	PullDate2=PullDate - pd.DateOffset(years=1)
	Mod16A2 = ee.ImageCollection("MODIS/006/MOD16A2");
	ModisShort = Mod16A2.filterDate(PullDate2, PullDate);
	ETShort = ModisShort.select('ET')
	ModisAnnualET = ETShort.sum()
	ModisAnnualET=ModisAnnualET.toDouble()
	ModisAnnualET=  ModisAnnualET.clip(AOI)
	PETShort = ModisShort.select('PET')
	ModisAnnualPET = PETShort.sum()
	ModisAnnualPET=ModisAnnualPET.toDouble()
	ModisAnnualPET=  ModisAnnualPET.clip(AOI)
	return(ModisAnnualET,ModisAnnualPET)
def RWS (MTBS_Row):
	shapefile =MTBS_Row
	features = []
	for i in range(shapefile.shape[0]):
		geom = shapefile.iloc[i:i+1,:] 
		jsonDict = eval(geom.to_json()) 
		geojsonDict = jsonDict['features'][0] 
		features.append(ee.Feature(geojsonDict)) 
	AOI = ee.FeatureCollection(features)
	StartDate=pd.to_datetime(str(int(MTBS_Row.Year))+'-'+str(int(MTBS_Row.StartMonth))+'-'+str(int(MTBS_Row.StartDay)))
	EndDate=StartDate + pd.DateOffset(days=3)
	Gridmet= ee.ImageCollection("IDAHO_EPSCOR/GRIDMET").filterDate(StartDate, EndDate)
	###Degrees Clockwise from north
	WindDir=Gridmet.select('th').mean().toDouble().clip(AOI)
	### Velocity at 10m
	WindSpd=Gridmet.select('vs').max().toDouble().clip(AOI)
	srtm = ee.Image('USGS/NED')
	slope = ee.Terrain.slope(srtm).clip(AOI)
	aspect=ee.Terrain.aspect(srtm).clip(AOI)
	Radians=slope.divide(ee.Number(180.0)).multiply(ee.Number(math.pi))
	RelativeWindDirection=WindDir.subtract(aspect).divide(ee.Number(180.0)).multiply(ee.Number(math.pi))
	combustionBuoyancy=10
	UaUb = WindSpd.divide(10) 
	### Set combustion bouyancy to 10 as it is unkwon 
	UaUb2=UaUb.multiply(UaUb)
	#2(UA/Ub)sinRadians*cos(slope)*Sin(relativeWindDirection)
	Secondterm=UaUb.multiply(2).multiply(Radians.sin()).multiply(RelativeWindDirection.cos())
	#Sin(slopeRadian)*2
	Thirdterm=Radians.sin().pow(2)
	EffWSPD=UaUb2.add(Secondterm.add(Thirdterm)).pow(.5)
	EffWSPD=EffWSPD.multiply(combustionBuoyancy).toDouble()
	EffWSPD_res=EffWSPD.reduceResolution(**{'reducer': ee.Reducer.mean()})
	return(EffWSPD,aspect,slope)
	
def GRIDMET(MTBS_Row,Variable,lagm=0):	
	TerraClimate=ee.ImageCollection("GRIDMET/DROUGHT")
	####
	shapefile =MTBS_Row
	features = []
	for i in range(shapefile.shape[0]):
		geom = shapefile.iloc[i:i+1,:] 
		jsonDict = eval(geom.to_json()) 
		geojsonDict = jsonDict['features'][0] 
		features.append(ee.Feature(geojsonDict)) 
	AOI = ee.FeatureCollection(features)

	if lagm>0:
		EndDate=pd.to_datetime(str(int(MTBS_Row.Year))+'-'+str(int(MTBS_Row.StartMonth))+'-'+str(int(MTBS_Row.StartDay)))
		StartDate=EndDate + pd.DateOffset(months=- lagm-1)
		Short = TerraClimate.filterDate(StartDate, EndDate).select(Variable).sum()
		OutputVar=Short.toDouble()
		OutputVar= OutputVar.clip(AOI)
	else:
		StartDate=pd.to_datetime(str(int(MTBS_Row.Year))+'-'+str(int(MTBS_Row.StartMonth))+'-'+str(int(MTBS_Row.StartDay)))
		#print(pd.DateOffset(month=5))
		EndDate=StartDate + pd.DateOffset(days=2)
		#print(StartDate,EndDate)
		Short = TerraClimate.filterDate(StartDate, EndDate).select(Variable).mean()
		OutputVar=Short.toDouble()
		OutputVar= Short.clip(AOI).toDouble()
	return(OutputVar)
	
def TerraClimate(MTBS_Row,Variable,lagm=0):	
	TerraClimate=ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")
	####
	shapefile =MTBS_Row
	features = []
	for i in range(shapefile.shape[0]):
		geom = shapefile.iloc[i:i+1,:] 
		jsonDict = eval(geom.to_json()) 
		geojsonDict = jsonDict['features'][0] 
		features.append(ee.Feature(geojsonDict)) 
	AOI = ee.FeatureCollection(features)

	if lagm>0:
		EndDate=pd.to_datetime(str(int(MTBS_Row.Year))+'-'+str(int(MTBS_Row.StartMonth))+'-'+str(int(MTBS_Row.StartDay)))
		StartDate=EndDate + pd.DateOffset(months=- lagm-1)
		Short = TerraClimate.filterDate(StartDate, EndDate).select(Variable).sum()
		OutputVar=Short.toDouble()
		OutputVar= OutputVar.clip(AOI)
	else:
		StartDate=pd.to_datetime(str(int(MTBS_Row.Year))+'-'+str(int(MTBS_Row.StartMonth))+'-'+str(int(MTBS_Row.StartDay)))
		#print(pd.DateOffset(month=5))
		EndDate=StartDate + pd.DateOffset(days=2)
		#print(StartDate,EndDate)
		Short = TerraClimate.filterDate(StartDate, EndDate).select(Variable).mean()
		OutputVar=Short.toDouble()
		OutputVar= Short.clip(AOI)
	return(OutputVar)	
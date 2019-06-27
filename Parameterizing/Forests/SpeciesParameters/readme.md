Species Parameterization
================

#### C:N ratios and Lignin

Most of the nitrogen, carbon and lignin parameters were started from existing LANDIS papers and the TRY database. The file NECN\_Folder contains the data taken from each paper and the database. Each Try data point is linked to an individual paper in the TRY files.

Some species were adapted from existing papers where the formula was assumed to have a carbon ratio of 1/2. This formula is, therefore, $C: N =.5 \\frac{g. nitrogen}{g. total}$. For species that were not included in either papers nor the plan database were researched individually, or given a value based on a genus/family/order similarity to a species in this or another Landis papers.

I started with values from previous LANDIS-II models and supplemented them with data from the TRY. This folder of accumulated values can be see [here](https://github.com/LANDIS-II-Foundation/Project-Southern-Appalachians-2018/blob/master/Parameterizing/Forests/SpeciesParameters/NECN_folder_12_20.xlsx).

Try data was queried by Species and by trait. In the even of multiple returns values were averaged.

Here is a example of some of the results obtained for Acer rubrum

|     |      sp     |                                     trait                                     |      unit     |   value   |
|-----|:-----------:|:-----------------------------------------------------------------------------:|:-------------:|:---------:|
| 5   | Acer rubrum |           Coarse woody debris (CWD) lignin content per CWD dry mass           |               | 77.750000 |
| 6   | Acer rubrum |                  Leaf nitrogen (N) content per leaf dry mass                  |               |     NA    |
| 7   | Acer rubrum | Branch coarse woody debris (CWD) nitrogen (N) content per branch CWD dry mass |       %       |  0.264200 |
| 8   | Acer rubrum |                  Leaf nitrogen (N) content per leaf dry mass                  |       %       |  1.752372 |
| 9   | Acer rubrum |                Stem (wood) carbon (C) content per stem dry mass               |       %       | 48.640000 |
| 10  | Acer rubrum |                              Crown (canopy) width                             |       cm      |  2.812360 |
| 11  | Acer rubrum |                        Leaf carbon/nitrogen (C/N) ratio                       | kg C per kg N | 24.700000 |
| 12  | Acer rubrum |                     Leaf carbon (C) content per leaf area                     |  kg C per m2  | 39.900000 |
| 13  | Acer rubrum |                  Leaf nitrogen (N) content per leaf dry mass                  |     mg g-1    | 17.359583 |
| 14  | Acer rubrum |                  Leaf nitrogen (N) content per leaf dry mass                  |    mg N g-1   | 19.623538 |

> TRY Data: Kattge, J., Bönisch, G., Günther, A., Wright, I., Zanne, A., Wirth, C., Reich, P.B. and the TRY Consortium (2012) TRY - Categorical Traits Dataset. Data from: TRY - a global database of plant traits. TRY File Archive <https://www.try-db.org/TryWeb/Data.php#3>

These values were used to update values for N:C and lignin, These records can be found at in the attached work. Values were also updated using records from:

> Davis, S. C., Dragan, K. E., Buyarski, C. R., & Thomas, R. B. (2009). High foliar and soil nitrogen concentrations in Central Appalachian forests. Ecosystems, 12(1), 46-56.

<img src="C:/Users/zjrobbin/Desktop/FS_meeting/FiguresMisc/Davis2009.PNG" width="60%" />

Values that were still missing were assessed as similar through genus then family and in one case by order. These decisions can be found in the document:

#### GGDmin,GDDmax,Frost,D3,FRT

These values were taken from the original linkages manual which can be found at <https://daac.ornl.gov/daacdata/model_archive/LINKAGES/comp/ORNL_TM-9519.pdf>. Other values were taken from existing LANDIS papers, these can be found in the NECN folder[here](https://github.com/LANDIS-II-Foundation/Project-Southern-Appalachians-2018/blob/master/Parameterizing/Forests/SpeciesParameters/NECN_folder_12_20.xlsx). Species that could not be found in either way were adapted from a qualitative assessment of range in comparison to known values for other species.

#### Functional Groups

In order to set up functional groups for the many species on the landscape, I wanted to compare their range of percipitation, termpature, elevation, and minimum vapor pressure. I did this using PRISM raster data on the 30 year normal. Using the imputation rasters from the Forest Service used in the [Initial Communities work](https://github.com/LANDIS-II-Foundation/Project-Southern-Appalachians-2018/tree/master/Parameterizing/Forests). I resampled the climate and elevation rasters to the same resolution and then, defined the minimum and maximum values for each species range (where it is &gt;.5 m2/ha) where the 20th and 80th quantile of each the total sample.

I then experimented with cluster groups by these features using a mean shift algorithm. Given a user-defined bandwidth, the algorithm finds the nth dimensional clustering of each group given a set of variables. I used different combinations of bandwidths and variables to find the clustering that seemed closest to a group of 3 functional groups. This is not always possible, if groups are too close or too far on nth dimensions, the natural cluster may be a different number. I then visually compared plots to see that these clusters made sense.

For variables, I decided on mean temperature, minimum precipitation, minimum vapor pressure deficit, and maximum elevation. Using mean shifting here is what the clusters looked like.

<img src="C:/Users/zjrobbin/Desktop/FS_meeting/Functional_Groups/Conifer1.jpg" width="40%" /> <img src="C:/Users/zjrobbin/Desktop/FS_meeting/Functional_Groups/Conifer2.jpg" width="40%" /> <img src="C:/Users/zjrobbin/Desktop/FS_meeting/Functional_Groups/Conifer3.jpg" width="40%" />

Additionally, I plotted each group on a PCA scale so that they could be visualized in one plot. Here the red lines represent the directions of the ordination and each conifer group is the colors on the right.

<img src="C:/Users/zjrobbin/Desktop/FS_meeting/Functional_Groups/ConiferPCA.jpg" width="40%" /> These three groups are labeled in future work.

Southern Pines: The Green group that has higher mean temp and higher min VPD. Northern Pines: The Red Group Colder temps, mid-elevation. Abies(Firs): The Blue group (one) High elevation and colder temps.

Here is the same thing for the hardwoods:

The hardwoods would not easily make 3 functional groups so I expanded it, this is because many of the species are niche specialists that either would have to be grouped as 2 big groups (which seemed coarse for so many species) or have many small individual groups.

<img src="C:/Users/zjrobbin/Desktop/FS_meeting/Functional_Groups/Hardwoods1 at1.2.jpg" width="40%" /> <img src="C:/Users/zjrobbin/Desktop/FS_meeting/Functional_Groups/Hardwoods2 at1.2.jpg" width="40%" /> <img src="C:/Users/zjrobbin/Desktop/FS_meeting/Functional_Groups/Hardwoods3 at1.2.jpg" width="40%" /> <img src="C:/Users/zjrobbin/Desktop/FS_meeting/Functional_Groups/HardwoodsPCA at1.2.jpg" width="40%" />

Here we see there are two major groups and then a lot of groups of outside groups. Looking to avoid parameterizing 8 functional groups, I took the visual grouping of group 5 with group 2 and clustering the right upper corner of the PCA together to create 4 groups:

This gave us the functional groups.:

Northern hardwoods

Southern Hardwoods

Riparians

Halesia

#### Max biomass

To find biomass and growth curves for parameterizing the NECN file of LANDIS-II for the Southern Apps Project. In LANDIS-II NECN, max biomass is a hypothetical maximum that a plot could hold if not in competition. To find that I am going to use known values of biomass and project them beyond a likely maximum.

This file is a continuation of the sorter that went through each FIA plot for the states of North Carolina, Tennessee, South Carolina, and Georgia, and calculated the Basal area for the plot and each species within it.

In this following loop, we will

Calculate the percent of total biomass Calculate the ratio of species biomass to total biomass We will isolate the top 60% stands in this ratio by age class, assuming this is growth under near to ideal conditions. Where this regression is at a 150% occupancy of a species, we will set as Max AGB Review graphs of this process and save the results

Bringing in the values from the FIA analysis. I isolated out the top 25% of plots by above ground carbon per age, assuming this to be the ideal growth for each species. These were then plotted as box plots which the comparison of each run could be simulated against, as a calibration measure.

#### Functional Group Parameters

Here we are parameterizing the growth curves of the LANDIS-II model against the box plots. Each line is a model run at a different point in the soil continuum from most sandy on the landscape to most clayey.

<img src="/Parameterizing/Forests/SpeciesParameters/QuerPrin/GrowthCurves.jpeg" width="40%" />

Using the FIA data we additionally fit a Mitscherlich curve to the relationship between age and biomass of the top 25 % of stands. This was compared to the biomass in an ideal scenario LANDIS-II runs.

<img src="/Parameterizing/Forests/SpeciesParameters/QuerPrin/QuerPrinAGB.jpeg" width="40%" />

We also compare the leaf area index of each species as a simulated forest, against established values for forests.

The values we are using for this are from:

> He, L., Chen, J. M., Pan, Y., Birdsey, R., & Kattge, J. (2012). Relationships between net primary productivity and forest stand age in US forests. Global Biogeochemical Cycles, 26(3).

On the left is a LANDIS output for chestnut oak and on the right are LAI values for a sampled oak/hickory forest.

<img src="/Parameterizing/Forests/SpeciesParameters/QuerPrin/LAI.jpeg" width="40%"/> 
<img src="/Parameterizing/Forests/SpeciesParameters/QuerPrin/Oak_Hickory.png" width="40%"/>

Here we looked at each of the major species in each functional group, as functional groups parameters in NECN are the major determinants of growth after the ideal growth ANPP parameter. If function group parameters seemed to work for the group but not an individual species, then ANPP was adjusted.

#### Northern Hardwoods

##### Quercus Prinus

    ## Warning: package 'gridExtra' was built under R version 3.4.4

(/QuerPrin/GrowthCurves.jpeg" width="40%")
<img src="/Parameterizing/Forests/SpeciesParameters/QuerPrin/QuerPrinAGB.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/QuerPrin/LAI.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/QuerPrin/Oak_Hickory.png" width="40%" />

##### Acer Rubra

<img src="/Parameterizing/Forests/SpeciesParameters/AcerRubr/GrowthCurves.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/AcerRubr/AcerRubrAGB.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/AcerRubr/LAI.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/AcerRubr/Maple_Beech_Birch.png" width="40%" />

#### Northern Conifers

##### Pinus Strobus

<img src="/Parameterizing/Forests/SpeciesParameters/PinuStorm/GrowthCurves.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/PinuStorm/PinuStormAGB.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/PinuStorm/LAI.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/PinuStorm/White_Pine.PNG" width="40%" />

##### Tsuga Canadensis

<img src="/Parameterizing/Forests/SpeciesParameters/Tsug Cand/GrowthCurves.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/Tsug Cand/Tsug CandAGB.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/Tsug Cand/LAI.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/Tsug Cand/Hemlock.PNG" width="40%" /> <img 
#### Southern Conifers

##### Pinus Virginiana

<img src="/Parameterizing/Forests/SpeciesParameters/PinuVirg/GrowthCurves.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/PinuVirg/PinuVirgAGB.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/PinuVirg/LAI.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/PinuVirg/Longleaf_Slash.png" width="40%" />

##### Pinus Taeda

<img src="/Parameterizing/Forests/SpeciesParameters/GrowthCurves.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/PinuTaedAGB.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/PinuTaed/LAI.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/Loblolly_Shortleaf.png" width="30%" />

#### Southern Hardwoods

##### liriodendron tulipifera

##### Carya Glabra

<img src="/Parameterizing/Forests/SpeciesParameters/caryglab/GrowthCurves.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/caryglab/CaryGlabAGB.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/caryglab/LAI.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/caryglab/Oak_Hickory.png" width="40%" />

#### Riparian

##### Betula Alleghaniensis

<img src="/Parameterizing/Forests/SpeciesParameters/BetuAlle/GrowthCurves.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/BetuAlle/BetuAlleAGB.jpeg" width="30%" /> <img src="/Parameterizing/Forests/SpeciesParameters/BetuAlle/LAI.jpeg" width="30%" /> <img src="/Parameterizing/Forests/SpeciesParameters/BetuAlle/Maple_Beech_Birch.png" width="30%" />

#### Abies

##### Frasier Fir

<img src="/Parameterizing/Forests/SpeciesParameters/Frasfir/GrowthCurves.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/Frasfir/FrasFirAGB.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/Frasfir/LAI.jpeg" width="40%" /> <img src="/Parameterizing/Forests/SpeciesParameters/Frasfir/Fir_Spruce.png" width="40%" />

If you are interested in the code to create this see the R-Markdown file on this page. To return to the main page click [here](https://github.com/LANDIS-II-Foundation/Project-Southern-Appalachians-2018)

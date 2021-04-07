Species Parameterization
================

#### C:N ratios and Lignin

Most of the nitrogen, carbon and lignin parameters were started from
existing LANDIS papers and the TRY database. The file NECN\_Folder
contains the data taken from each paper and the database. Each Try data
point is linked to an individual paper in the TRY files.

Some species were adapted from existing papers where the formula was
assumed to have a carbon ratio of 1/2. This formula is, therefore,
\(C: N =.5 \frac{g. nitrogen}{g. total}\). For species that were not
included in either papers nor the plan database were researched
individually, or given a value based on a genus/family/order similarity
to a species in this or another Landis papers.

I started with values from previous LANDIS-II models and supplemented
them with data from the TRY. This folder of accumulated values can be
see
[here](https://github.com/LANDIS-II-Foundation/Project-Southern-Appalachians-2018/blob/master/Parameterizing/Forests/SpeciesParameters/NECN_folder_12_20.xlsx).

Try data was queried by Species and by trait. In the even of multiple
returns values were averaged.

Here is a example of some of the results obtained for Acer rubrum

|     |     sp      |                                     trait                                     |     unit      |   value   |
| :-- | :---------: | :---------------------------------------------------------------------------: | :-----------: | :-------: |
| 1   | Acer rubrum |           Coarse woody debris (CWD) lignin content per CWD dry mass           |               | 77.750000 |
| 9   | Acer rubrum |                  Leaf nitrogen (N) content per leaf dry mass                  |               |    NA     |
| 16  | Acer rubrum | Branch coarse woody debris (CWD) nitrogen (N) content per branch CWD dry mass |       %       | 0.264200  |
| 23  | Acer rubrum |                  Leaf nitrogen (N) content per leaf dry mass                  |       %       | 1.752372  |
| 59  | Acer rubrum |               Stem (wood) carbon (C) content per stem dry mass                |       %       | 48.640000 |
| 73  | Acer rubrum |                             Crown (canopy) width                              |      cm       | 2.812360  |
| 98  | Acer rubrum |                       Leaf carbon/nitrogen (C/N) ratio                        | kg C per kg N | 24.700000 |
| 111 | Acer rubrum |                     Leaf carbon (C) content per leaf area                     |  kg C per m2  | 39.900000 |
| 137 | Acer rubrum |                  Leaf nitrogen (N) content per leaf dry mass                  |    mg g-1     | 17.359583 |
| 158 | Acer rubrum |                  Leaf nitrogen (N) content per leaf dry mass                  |   mg N g-1    | 19.623538 |

> TRY Data: Kattge, J., B?nisch, G., G?nther, A., Wright, I., Zanne, A.,
> Wirth, C., Reich, P.B. and the TRY Consortium (2012) TRY - Categorical
> Traits Dataset. Data from: TRY - a global database of plant traits.
> TRY File Archive <https://www.try-db.org/TryWeb/Data.php#3>

These values were used to update values for N:C and lignin, These
records can be found at in the attached work. Values were also updated
using records from:

> Davis, S. C., Dragan, K. E., Buyarski, C. R., & Thomas, R. B. (2009).
> High foliar and soil nitrogen concentrations in Central Appalachian
> forests. Ecosystems, 12(1), 46-56.

Values that were still missing were assessed as similar through genus
then family and in one case by order. These decisions can be found in
the document
[here](https://github.com/LANDIS-II-Foundation/Project-Southern-Appalachians-2018/blob/master/Parameterizing/Forests/SpeciesParameters/Species%20by%20Species%20parameters_Documentation.xlsx)

#### GGDmin,GDDmax,Frost,D3,FRT

These values were taken from the original linkages manual which can be
found at
<https://daac.ornl.gov/daacdata/model_archive/LINKAGES/comp/ORNL_TM-9519.pdf>.
Other values were taken from existing LANDIS papers, these can be found
in the NECN
folder[here](https://github.com/LANDIS-II-Foundation/Project-Southern-Appalachians-2018/blob/master/Parameterizing/Forests/SpeciesParameters/NECN_folder_12_20.xlsx).
Species that could not be found in either way were adapted from a
qualitative assessment of range in comparison to known values for other
species.

#### Functional Groups

In order to set up functional groups for the many species on the
landscape, I wanted to compare their range of percipitation, termpature,
elevation, and minimum vapor pressure. I did this using PRISM raster
data on the 30 year normal. Using the imputation rasters from the Forest
Service used in the [Initial Communities
work](https://github.com/LANDIS-II-Foundation/Project-Southern-Appalachians-2018/tree/master/Parameterizing/Forests).
I resampled the climate and elevation rasters to the same resolution and
then, defined the minimum and maximum values for each species range
(where it is \>.5 m2/ha) where the 20th and 80th quantile of each the
total sample.

I then experimented with cluster groups by these features using a mean
shift algorithm. Given a user-defined bandwidth, the algorithm finds the
nth dimensional clustering of each group given a set of variables. I
used different combinations of bandwidths and variables to find the
clustering that seemed closest to a group of 3 functional groups. This
is not always possible, if groups are too close or too far on nth
dimensions, the natural cluster may be a different number. I then
visually compared plots to see that these clusters made sense.

For variables, I decided on mean temperature, minimum precipitation,
minimum vapor pressure deficit, and maximum elevation. Using mean
shifting here is what the clusters looked like.

These three groups are labeled in future work.

Southern Pines: The Green group that has higher mean temp and higher min
VPD. Northern Pines: The Red Group Colder temps, mid-elevation.
Abies(Firs): The Blue group (one) High elevation and colder temps.

We did the same thing for the hardwoods:

The hardwoods would not easily make 3 functional groups so I expanded
it, this is because many of the species are niche specialists that
either would have to be grouped as 2 big groups (which seemed coarse for
so many species) or have many small individual groups.

Looking to avoid parameterizing 8 functional groups, I took the visual
grouping of 3 groups and clustering the right upper corner of the PCA
together to create 3 groups:

This gave us the functional groups:

Northern hardwoods

Southern Hardwoods

Riparians

#### Max biomass

To find biomass and growth curves for parameterizing the NECN file of
LANDIS-II for the Southern Apps Project. In LANDIS-II NECN, max biomass
is a hypothetical maximum that a plot could hold if not in competition.
To find that I am going to use known values of biomass and project them
beyond a likely maximum.

This file is a continuation of the sorter that went through each FIA
plot for the states of North Carolina, Tennessee, South Carolina, and
Georgia, and calculated the Basal area for the plot and each species
within it.

In this following loop, we will

Calculate the percent of total biomass Calculate the ratio of species
biomass to total biomass We will isolate the top 60% stands in this
ratio by age class, assuming this is growth under near to ideal
conditions. Where this regression is at a 150% occupancy of a species,
we will set as Max AGB Review graphs of this process and save the
results

Bringing in the values from the FIA analysis. I isolated out the top 25%
of plots by above ground carbon per age, assuming this to be the ideal
growth for each species. These were then plotted as box plots which the
comparison of each run could be simulated against, as a calibration
measure.

## Functional Groups

To fit the parameters around the functional group we compare the most
prevalent species to the FIA record of biomass per site, and a published
value for LAI. For each, we ran single-cell simulations using the mean
landscape values for NECN inputs. Each stand was initialized with two
cohorts (one age 0 and one age 10) both with densities of 1000 g/m2.

#### Northern Conifer

For northern conifers the most prevelent species was Pinus strobus.

![Northern Conifer
Biomass](Functional%20Group%20Validating/Northern%20Conifer%20Biomass.png)

![Northern Conifer
Biomass](Functional%20Group%20Validating/Northern%20Conifer%20LAI.png)

Compared to white/red/jack pine values from He, L., Chen, J. M., Pan,
Y., Birdsey, R., & Kattge, J. (2012). Relationships between net primary
productivity and forest stand age in US forests. Global Biogeochemical
Cycles, 26(3).

![White Pine
LAI](Functional%20Group%20Validating/PicsofTrees/White_Pine.png)

### Southern Conifer

For southern conifers we compared published values from Pinus Taeda

![Pinus taeda
Biomass](Functional%20Group%20Validating/Southern%20Conifer%20Biomass.png)

![Pinus Taeda
LAI](Functional%20Group%20Validating/Southern%20Conifer%20LAI.png)
Here are values from Loblolly/shortleaf pine stands from He, L., Chen,
J. M., Pan, Y., Birdsey, R., & Kattge, J. (2012). Relationships between
net primary productivity and forest stand age in US forests. Global
Biogeochemical Cycles, 26(3).

![Lobolloy
LAI](Functional%20Group%20Validating/PicsofTrees/Loblloly_Shortleaf.png)

\#\#\#Southern Hardwoods For southern hardwoods we compared Tulip poplar
(Liriodendron tulipifera).

![Tulip poplar
Biomass](Functional%20Group%20Validating/Southern%20Hardwoods%20Biomass.png)
![Tulip poplar
LAI](Functional%20Group%20Validating/Southern%20Hardwoods%20LAI.png)

Published values from Mixed Tulip-poplar/ Sweet gum forests

Brown, M. J., & Parker, G. G. (1994). Canopy light transmittance in a
chronosequence of mixed-species deciduous forests. Canadian Journal of
Forest Research, 24(8), 1694-1703.

![Tupli poplar
LAI](Functional%20Group%20Validating/PicsofTrees/Tulip_poplar.jpg)

### Northern Hardwoods

Chestnut oak (Quercus montana) was used to compare northern hardwoods

![Chestnut Oak
Biomass](Functional%20Group%20Validating/Northern%20Hardwood%20Biomass.png)

![Chestnut Oak
Biomass](Functional%20Group%20Validating/Northern%20Hardwood%20LAI.png)

Published Chestnut Oak LAI He, L., Chen, J. M., Pan, Y., Birdsey, R., &
Kattge, J. (2012). Relationships between net primary productivity and
forest stand age in US forests. Global Biogeochemical Cycles, 26(3).

![Chestnut Oak
LAI](Functional%20Group%20Validating/PicsofTrees/Oak_Hickory.png)

### Hemlock

Hemlock is its own functional group.

![Hemlock
Biomass](Functional%20Group%20Validating/Hemlock%20Biomass.png)

![Hemlock
LAI](Functional%20Group%20Validating/Hemlock%20LAI.png)

![Hemlock
LAI](Functional%20Group%20Validating/PicsofTrees/Hemlock.png)

### Swampy

For this functional group we used Betula Alleghaniensis

![Yellow
Birch](Functional%20Group%20Validating/Swampy_Biomass.png)

![Yellow
Birch](Functional%20Group%20Validating/Swampy_LAI.png)

We compared this to Maple, Beech and Birch stands

![Beech
LAI](Functional%20Group%20Validating/PicsofTrees/Maple_Beech_Birch.png)

### Abies

For this functional group we used Fraser Fir, unfortunately there were
very few FIA plots to use as comparison, so we used the max biomass
esablished in prior LANDIS-II studies of the species.

![Fraser
fir](Functional%20Group%20Validating/Abies_Biomass.png)

![Fraser
fir](Functional%20Group%20Validating/Abies_LAI.png)

Published Fir-Spruce stand LAI He, L., Chen, J. M., Pan, Y., Birdsey,
R., & Kattge, J. (2012). Relationships between net primary productivity
and forest stand age in US forests. Global Biogeochemical Cycles, 26(3).

![FF
LAI](Functional%20Group%20Validating/PicsofTrees/Fir_Spruce.png)

If you are interested in the code to create this see the R-Markdown file
on this page. To return to the main page click
[here](https://github.com/LANDIS-II-Foundation/Project-Southern-Appalachians-2018)

Rx\_Burn\_Data
================
Kate Jones
updated 01/21/2021

This markdown describes acquiring, cleaning, and processing RX data in
Western North Carolina to create an RX input raster for the S. Apps
SCRPPLE extension.

**Data**

The data to create the RX surface were predominantly gathered through
personal correspondence with managers in western North Carolina. We were
unable to acquire data for TN, GA, or SC, so the data values calculated
for NC were extrapolated to the other states in our study area. Data to
estimate federal forest burning was gathered separately from a USFS
report - no spatial data for Federal lands were acquired, so we expect
RX burning on federal lands is underestimated in this work flow. It
should also be noted, that the numbers reflected in the USFS report are
annual targets, which may or may not be met each year, so we handled
these data separately than the observed RX data gathered from the NC
agencies.

Contacts where RX data were acquired:

**1) NC Forest Service (NCFS)**  
- Christian D. Vose \| Information Technology Branch Head \|
<christian.vose@ncagr.gov>

**2) NC State Parks (NCSP)**  
- John Amoroso \| Planning & GIS \| <john.amoroso@ncparks.gov>

**3) NC Wildlife Resources Commission (NCWRC)**  
- Ryan Jacobs \| Wildlife Forest Manager \| <Ryan.Jacobs@ncwildlife.org>

**Step 1: Load & Project Data**

The data collected from each agency were in slightly different formats,
but all spatial. We read all the datasets in and project to UTM 17N and
crop them to the AOI.

**Step 2: Plot and Check**

To make sure everything spatially aligns and looks as expected, we plot
the data. Looks ok, but we’ll want to zoom in spatially to the AOI and
change plotting for the polygons layers.

\*The NCWRC and NC State Parks provided area data, while the NCFS
provided point data. We will convert all data to points (using area
centroids) for consistency.

![](Rx_data_processing_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

**Step 3: Better Plotting**

Once zoomed into the AOI, we can better see the RX locations for the
three agencies. You’ll notice there are hard cutoffs in the AOI between
areas with ignitions and areas without ignitions - this is because we
were only able to collect data from NC, so these boundaries follow state
lines. We can see the RX burns are fairly clustered (i.e. the
observations generally occur in the same locations, which are contained
within management boundaries). The exception is the NCFS in grey, which
had an unexpected number of burns that appeared on private lands
(according to the land ownership data we use). This may explain why the
NCFS RX ignitions are more distributed throughout the state and not
confined to certain land ownership. More on this later.

![](Rx_data_processing_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

**Step 4: Dividing Data by Jurisdiction/Ownership**

There could be multiple iterations of this process, but for now, we’ve
grouped the land ownership layer into Federal, State, Private, and Local
(local is likely an unnecessary class, but without further discussion on
where to group it, it stands alone). The land ownership data will be
intersected with the RX locations.

Land ownership data was acquired via correspondence with NIFCC and
downloaded from this link:
<https://wfdss.usgs.gov/wfdss_help/WFDSSHelp_SMA_Layer.html>, then
navigate to the Data Downloads page.

![](Rx_data_processing_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

**Step 5: Format Rx Datasets as Points** For now, converting everything
to points (i.e. centroids of polygons) so it’s in the same format and
assigning the polygon dataframe (attributes) to the centroids.

**Step 6: Getting all of the data in the same format**

Across the datasets, there are different formats for dates, burn names,
objectids, etc. Need to get everything into a format that has the
following columns:

-OBJECTID  
-Burn Unit (if possible)  
-Year  
-Month  
-Date  
-Acreage

To get the above columns in a unified format, the WRC data needs the
most cleaning. Cleaning steps commented in code.

**Step 7: Final data formatting**

Combine all formatted data into single spatial object with tabular
information.

**Quick Diagnostic: Monthly RX occurrences all years, all data sources**

    ## [1]  1 12

![](Rx_data_processing_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

**Quick Diagnostic: Historic RX Acreage Burned by Ownership**

![](Rx_data_processing_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

**Quick Diagnostic: Historic Burn Sizes by Ownership**

![](Rx_data_processing_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->![](Rx_data_processing_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

**Step 8: Determine Burned Area %’s by Land Ownership**

**Quick Diagnostic: Annual Acres Burned by Ownership**

![](Rx_data_processing_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

![](Rx_data_processing_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

**Step 9: Assigning raster cells a percentage of total acres burned
across the landscape**

This represents the probability that a given cell will burn. The tables
below each raster show the relative burn probabilities and the
hypothetical burned areas used to calculate the probabilities for each
ownership for the 10/50/100 scenarios.

![](Rx_data_processing_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

| fspl    | total\_burnedacres\_injuris | total\_burnedacres\_landscape | perc\_of\_land\_burn | round\_perc\_land\_burn |
|:--------|----------------------------:|------------------------------:|---------------------:|------------------------:|
| Fed     |                       80000 |                        141153 |            0.5667609 |                    0.57 |
| Local   |                         126 |                        141153 |            0.0008926 |                    0.00 |
| Private |                       18372 |                        141153 |            0.1301566 |                    0.13 |
| State   |                       42655 |                        141153 |            0.3021898 |                    0.30 |

This table shows the values (empirical and derived) to represent RX fire
across the landscape when the NF is meeting 10% of their acreage goals.
All RX values are held constant across the 10/50/100 % scenarios, as
these are observed data.

![](Rx_data_processing_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

| fspl    | total\_burnedacres\_injuris | total\_burnedacres\_landscape | perc\_of\_land\_burn | round\_perc\_land\_burn |
|:--------|----------------------------:|------------------------------:|---------------------:|------------------------:|
| Fed     |                      400000 |                        461153 |            0.8673911 |                    0.87 |
| Local   |                         126 |                        461153 |            0.0002732 |                    0.00 |
| Private |                       18372 |                        461153 |            0.0398393 |                    0.04 |
| State   |                       42655 |                        461153 |            0.0924964 |                    0.09 |

This table shows the values (empirical and derived) to represent RX fire
across the landscape when the NF is meeting 50% of their acreage goals.
All RX values are held constant across the 10/50/100 % scenarios, as
these are observed data.

![](Rx_data_processing_files/figure-gfm/unnamed-chunk-16-3.png)<!-- -->

| fspl    | total\_burnedacres\_injuris | total\_burnedacres\_landscape | perc\_of\_land\_burn | round\_perc\_land\_burn |
|:--------|----------------------------:|------------------------------:|---------------------:|------------------------:|
| Fed     |                      800000 |                        861153 |            0.9289871 |                    0.93 |
| Local   |                         126 |                        861153 |            0.0001463 |                    0.00 |
| Private |                       18372 |                        861153 |            0.0213342 |                    0.02 |
| State   |                       42655 |                        861153 |            0.0495324 |                    0.05 |

This table shows the values (empirical and derived) to represent RX fire
across the landscape when the NF is meeting 100% of their acreage goals.
All RX values are held constant across the 10/50/100 % scenarios, as
these are observed data.

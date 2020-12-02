# megaSDM
## Description
This package can efficiently create and project species distribution models using the MaxEnt framework and parallel processing. It can find and download occurrence data for a list of species on GBIF (Global Biodiversity Information Facility), environmentally subsample the occurrences to mitigate spatial bias, generate background (pseudo-absence) points, train the model and project it to different times (incorporating dispersal rate of each species and intermediate range fluctuations), and create species richness maps for each time period and taxon. 

## Contributors

Benjamin Shipley

Renee Bach

Younje Do

Heather Strathearn

Jenny McGuire

Bistra Dilkina

## Dependencies *(updated 18 November 2020)*
#### maxent.jar file (may be downloaded at https://github.com/mrmaxent/Maxent)
#### R version 3.6 or greater (last tested with v. 4.0.2)
### R Package (Version Number of Last Tested)
dplyr	(1.0.2)

gtools	(3.8.2)

plotfunctions	(1.4)

raster	(3.3.13)

rgbif	(3.3.0)

rgdal	(1.5.18)

rgeos	(0.5.5)

sampSurf	(0.7.5)

sp (1.4.2)

## Installation Instructions
Refer to the provided worked example (https://github.com/brshipley/megaSDM/blob/master/EXAMPLE/WorkedExample.R) for detailed instructions on the setup and configuration of this program.
## Usage Example
We have provided an example (https://github.com/brshipley/megaSDM/blob/master/EXAMPLE.zip) using occurrence data from GBIF for 6 North American mammal species. Environmental data are from the WorldClim database (Hijmans et al. 2005; https://www.worldclim.org/), and data on dispersal rate (in km/year) were collected by HS. The R script "WorkedExample.R" displays the entire functionality of megaSDM in a cohesive workflow, from data collection to the analysis and presentation of results. 

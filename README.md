# megaSDM
## Description
This package can efficiently create and project species distribution models using the MaxEnt framework and parallel processing. It can find and download occurrence data for a list of species on GBIF (Global Biodiversity Information Facility), environmentally subsample the occurrences to mitigate spatial bias, generate background (pseudo-absence) points, train the model and project it to different times (incorporating dispersal rate of each species and intermediate range fluctuations), and create species richness maps for each time period and taxon. 

## Contributors

Benjamin R. Shipley 
Renee Bach
Daniel Do
Heather Strathearn
Bistra Dilkina
Jenny L. McGuire

## Dependencies *(updated 16 Mar 2020)*
#### maxent.jar file (may be downloaded at https://github.com/mrmaxent/Maxent)
#### R version 3.6.1
### R Package (Version Number)

dplyr	(0.8.1)

gtools	(3.8.1)

parallel	(3.5.1)

plotfunctions	(1.3)

raster	(2.9-5)

rgbif	(1.3)

rgdal	(1.4-3)

rgeos	(0.4-3)

sampSurf	(0.7-4)

## Installation Instructions
Refer to the provided setup guide (https://github.com/brshipley/megaSDM/blob/master/megaSDM_Setup.pdf) for detailed instructions on the setup and configuration of this program.
## Usage Example
We have provided an example (https://github.com/brshipley/megaSDM/blob/master/EXAMPLE.zip) using occurrence data from GBIF for 6 North American mammal species. Environmental data are from the WorldClim database (Hijmans et al. 2005; https://www.worldclim.org/version1), and data on dispersal rate (in km/year) were collected by HS. Note that the specific file paths described in the configuration file within the .zip folder may need to be modified once downloaded.

## Stand-Alone Functions
In addition to the main workflow, this repository includes several stand-alone functions that reproduce some aspects of the primary scripts without the need for the configuration file or prior file management. A README file describing each function may be found at (https://github.com/brshipley/megaSDM/blob/master/Stand-Alone/README.md), along with the scripts of the functions.

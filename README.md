# megaSDM
## Description
This package can efficiently create and project species distribution models using the MaxEnt framework and parallel processing. It can find and download occurrence data for a list of species on GBIF (Global Biodiversity Information Facility), environmentally subsample the occurrences to mitigate spatial bias, generate background (pseudo-absence) points, train the model and project it to different times (incorporating dispersal rate of each species and intermediate range fluctuations), and create species richness maps for each time period and taxon. 

##Contributors

## Dependencies *(updated 9 Sept 2019)*
#### maxent.jar file (may be downloaded at https://github.com/mrmaxent/Maxent)
#### R version 3.6.1
### R Package (Version Number)
biomod2	(3.3-7.1)

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
Refer to the provided setup guide (https://github.com/brshipley/megaSDM/megaSDM_Setup.pdf) for detailed instructions on the setup and configuration of this program.
## Usage Example
We have provided an example (____) using occurrence data from GBIF for 6 North American vertebrate species (3 reptiles, 3 mammals) and selected environmental data from the WorldClim database (Hijmans et al. 2005; https://www.worldclim.org/version1). The configuration file used in this example may be found at (_____). Note that the specific file paths described in the configuration file may need to be modified once downloaded.

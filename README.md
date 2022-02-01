
<!-- README.md is generated from README.Rmd. Please edit that file -->

# megaSDM

<!-- badges: start -->

<!-- badges: end -->

## Description

This package can efficiently create and project species distribution
models using the MaxEnt framework and parallel processing. It can find
and download occurrence data for a list of species on GBIF (Global
Biodiversity Information Facility), environmentally subsample the
occurrences to mitigate spatial bias, generate background
(pseudo-absence) points, train the model and project it to different
times (incorporating dispersal rate of each species and intermediate
range fluctuations), and create species richness maps for each time
period and taxon.

## Contributors

Benjamin Shipley

Renee Bach

Younje Do

Heather Strathearn

Jenny McGuire

Bistra Dilkina

## Dependencies *(updated 31 January 2022)*

#### maxent.jar file (may be downloaded at <https://github.com/mrmaxent/Maxent>)

#### R version 3.6 or greater (last tested with v. 4.1.2)

### R Package (Version Number of Last Tested)

dplyr (1.0.7)

gtools (3.9.2)

plotfunctions (1.4)

raster (3.5.15)

rgbif (3.6.0)

rgdal (1.5.28)

rgeos (0.5.9)

sampSurf (0.7.6)

sp (1.4.6)

## Vignette

We have provided an example vignette
(brshipley/megaSDM/megaSDM\_vignette.html) using occurrence data from
GBIF for 6 North American mammal species. Environmental data are from
the WorldClim database (Hijmans et al. 2005;
<https://www.worldclim.org/>), and data on dispersal rate (in km/year)
were collected by HS. The html file (and acossiated R Markdown file
within the package) displays the entire functionality of megaSDM in a
cohesive workflow, from data collection to the analysis and presentation
of results.

## Installation

In R, use `devtools::install_github("brshipley/megaSDM", build_vignettes
= TRUE)` to install the package with the vignette (see above) included.
To access the vignette itself, run `??megaSDM_vignette`, and the
vignette will be loaded in the HTML viewer.

## Reference and Citation

For more information about the methods employed in the package and highlighted features, 
refer to <https://doi.org/10.1111/ecog.05450>. 

Please cite this package as:
Shipley, B. R., Bach, R., Do, Y., Strathearn, H., McGuire, J. L., & Dilkina, B. (2022).
megaSDM: integrating dispersal and time‐step analyses into species distribution models.
Ecography, 2022: e05450.


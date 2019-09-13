# Stand-Alone Functions

## EnvCorrelation
### Description:
A wrapper for raster::layerStats (Hijmans et al. 2017) that returns a correlation matrix for a set of rasters and warning messages if two climate layers are highly correlated.
### Parameters:
*envlayers* (list or RasterStack): A list or RasterStack of rasters to be correlated

*threshold* (numeric): if correlation coefficients are above the given value, a warning message is printed to the console.

## VarelaSample
### Description:
Environmentally subsamples given occurrence points based on a set of climate rasters (see Varela et al. 2014, Castellanos et al. 2018) and returns either a data frame of subsampled occurrences or a data frame showing the number of subsampled occurrences generated for each given number of bins. Climate binning can be done on the climate variables themselves or on a specified number of principal component axes. 
### Parameters:
*occurrences* (data.frame or table): Species Occurrences with columns of coordinates [long, lat] in decimal degree format

*env* (RasterLayer or Stack): A RasterLayer or RasterStack object of all environmental data

*nbins* (numeric) : number of bins to divide each environmental layer or PC axis into for stratified subsampling. If nbins is a range of values, a data frame of the number of subsampled points for the range of bin numbers will be output instead of the subsampled occurrences

*PCA* (logical (Y/N)): . Should environmental data be binned by Principal Component Axes instead of the actual data values (Castellanos et al. 2018)?

*PCAxes* (numeric): The number of PC Axes to use for the environmental subsampling. if not given, the number of axes required to explain 95% of the data variance will be used.

## BackgroundPoints
### Description:
Generates and (optionally) environmentally subsamples a set of background points (see VarelaSample). If a shapefile is included, 1/2 of the background points will be generated from within the shapefile, providing a spatial weighting to the background points (see Fig. 5 in this paper)
### Parameters:
*nbg* (numeric) : The number of background points wanted (if method = “Varela”, this number will not be exact)

*envdata* (RasterLayer or Stack): A RasterLayer or RasterStack object of all environmental data

*buffer* (shapefile): a shapefile polygon (usually a buffer around the occurrence points) for spatial weighting

*method* (string: “Varela” or “random”): How should background points be generated. If method = “Varela”, environmental subsampling will occur.

*nbins* (numeric): number of bins to divide each environmental layer or PC axis into for stratified subsampling.

*PCA* (numeric or logical (Y/N)): Should environmental data be binned by Principal Component Axes instead of the actual data values (Castellanos et al. 2018)? If a number [n] is given, n PC Axis will be used for subsampling. Otherwise, the number of axes required to explain 95% of the data variance will be used.

## dispersalRateExp
### Description:
Given a mean or median dispersal rate per unit time and the number of time units desired, calculates the mean and median total distance dispersed and (optionally) generates a graph of the dispersal distribution and a vector of probability values.
### Parameters:
*dist* (numeric): the distance dispersed in one time step (years,generations,etc.)

*time_steps* (numeric): the number of time steps to model dispersal through (>0)

*method* (string: "mean" or "median"): Is the given dispersal rate the mean or a median rate?

*value* (numeric vector): A vector of distances. If given, probability of dispersal to these distances will be calculated

*plot*: (logical). Should a plot of the total dispersal probability be generated?

## RichnessMaps
### Description:
Generates consistent species richness maps for binary distribution rasters
### Parameters:
*speciesrasters* (list or RasterStack): A list or RasterStack of binary rasters (1 for presence, 0 for absence). 

*plotmap* (logical): Should a species richness map be displayed in R?

*pdffile* (string): a character string of a file name for a PDF output (must include ".pdf")

*title* (string): a character string of the PDF and/or displayed map title

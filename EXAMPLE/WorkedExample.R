#### Worked example for megaSDM ####
#### This example demonstrates much of the functionality of megaSDM
#### and can be used as a template for other SDM analyses.

# Setup---------------------------
# the modelling portion fo megaSDM relies on the MaxEnt algorithm,
# so the "maxent.jar" must be downloaded before modelling can occur.

# The "maxent.jar" file can be found at https://github.com/mrmaxent/Maxent

# In this example, the "maxent.jar" file is located in "./EXAMPLE/TestRun/models"

# Double-click the "maxent.jar" file to ensure that it can open.
# If the maxent.jar dialog box does not load when opened (this is common when using MacOS),
# the Java Runtime Environment may need to be downloaded before running megaSDM
# (https://www.oracle.com/java/technologies/javase-jre8-downloads.html).

# File Management:

# Copy the folder with the data and the WorkedExample into a new directory
# (e.g., I used the desktop folder "C:/Users/bshipley6/Desktop/EXAMPLE")

# Next, copy the folder with the megaSDM package in it into a directory
# (e.g., "C:/Users/bshipley6/Desktop/EXAMPLE/megaSDM")
  #NOTE keep the name of this folder the same so that installation of megaSDM can proceed


#Downloading the package:
#Install and Load these two packages:
install.packages(c("devtools", "roxygen2"))
library(devtools)
library(roxygen2)

#Set the working directory to the package folder and install
setwd("C:/Users/bshipley6/Desktop/EXAMPLE")
install("megaSDM")

#When R asks if packages should be updated, enter a blank line.
library(megaSDM)

# 1. Training/Study Env Projection/Clipping-----------------------------

# Species distribution models use a training area and a study area. The training area is the area where
# the model will be trained on (i.e., where the occurrence and background points for the model generation
# are located). The study area is the region of interest (i.e., where the model will be projected and habitat
# suitability will be predicted for both the current time period and future/past time periods).

# The first function "TrainStudyEnv" creates and manages the environmental rasters of the current time period for
# the desired training and study areas.

# First, define "envoutput", which is where the training and study rasters will be printed out to.
# In this example it will be the same location the TestRun will be located in.
envoutput <- "C:/Users/bshipley6/Desktop/EXAMPLE/TestRun"

# Make a list of environmental rasters that we will use to clip, resample, and reproject
# In this example, navigate to the current environmental files within "Data" (/trainingarea)
input_TA <- list.files("C:/Users/bshipley6/Desktop/EXAMPLE/Data/trainingarea", pattern = ".bil$", full.names = TRUE)

# For this example, our training and study areas will be exactly the same, and we don't need
# to reproject or resample
TSEnv <- TrainStudyEnv(input_TA = input_TA,
                       output = envoutput,
                       clipTrain = c(-91.5, -75, 25.5, 36),
                       clipStudy = c(-91.5, -75, 25.5, 36))
# By saving this as an object "TSEnv", the raster are saved within the R environment for use later

# 2. Forecasted/Hindcasted Env Projections-------------------------------

# After getting the training and study environmental layers for the current time period,
# the next function ("PredictEnv") clips, resamples, and reprojects the Future/Past environmental
# layers to the characteristics of the given study region.

# This function requires lists of the forecasted/hindcasted environmental files.
# Navigate to the directory with RCP4.5 and click on a file in the correct time period for each
Env2050_4.5 <- list.files(dirname(file.choose()), pattern = ".bil$", full.names = TRUE)
Env2070_4.5 <- list.files(dirname(file.choose()), pattern = ".bil$", full.names = TRUE)
Env4.5 <- list(Env2050_4.5, Env2070_4.5)

# This function takes the study area rasters of the current time period and clips and projects
# the forecasted/hindcasted environmental rasters to match.

# The "time_periods" argument contains the current time (the time of the training and study rasters)
#  first and then the time periods for the forecast/hindcast

PredictEnv(studylayers = TSEnv$study,
           futurelayers = Env4.5,
           time_periods = c(2010, 2050, 2070),
           output = envoutput,
           scenario_name = "RCP4.5")

# Repeat for Scenario RCP8.5 (navigate there)
Env2050_8.5 <- list.files(dirname(file.choose()), pattern = ".bil$", full.names = TRUE)
Env2070_8.5 <- list.files(dirname(file.choose()), pattern = ".bil$", full.names = TRUE)

Env8.5 <- list(Env2050_8.5, Env2070_8.5)
PredictEnv(studylayers = TSEnv$study,
           futurelayers = Env8.5,
           time_periods = c(2010, 2050, 2070),
           output = envoutput,
           scenario_name = "RCP8.5")

# 3. Occurrence Collection--------------------------
# This function gathers occurrence points from GBIF, given a list of species and an extent

# The extent should be the same (or similar to) as the extent of the training area.
# Given in latitude/longitude coordinates:
extent_occ <- c(-91.5, -75, 25.5,36)

# A list of southeastern mammals for this example
spplist <- c("Puma concolor coryi",
             "Podomys floridanus",
             "Sylvilagus aquaticus",
             "Sylvilagus palustris",
             "Geomys pinetis",
             "Neofiber alleni")

# Define the file folder where the occurrences will be written
# (if this folder doesn't already exist, megaSDM will make it)
occ_output <- "C:/Users/bshipley6/Desktop/EXAMPLE/TestRun/occurrences"

# This function only takes occurrences from the described trainingarea extent.
Occurrences <- OccurrenceCollection(spplist = spplist,
                                    output = occ_output,
                                    trainingarea = extent_occ)

# Because one species was renamed, rename species list to reflect taxonomy changes
spplist <- Occurrences$Scientific.Name

# 4. Occurrence Management--------------------------
# After collecting the occurrences, they must be formatted so that they are consistent and able
# to be read within MaxEnt. This function takes a list of occurrence files, extracts environmental
# data at each point, and if necessary, subsamples the dataset for more accurate modelling.

# First, get the list of the occurrence files
occlist <- list.files(occ_output, pattern = ".csv", full.names = TRUE)

# If the output is set to the same place as in the previous function, it will overwrite the occurrence files
OccurrenceManagement(occlist = occlist,
                     output = occ_output,
                     envextract = TRUE,
                     envsample = TRUE,
                     nbins = 25,
                     envdata = TSEnv$training)

# 5. Background Buffers-----------------------------

# Spatially-constrained background poitns can be more effective for SDMs. To constrain background points,
# megaSDM uses buffers to sample points within a certain radius of the occurrence points. This function
# generates those buffers.

# We need to get the list of occurrence files again (even if they were written out in the same folder as before)
occlist <- list.files(occ_output, pattern = ".csv", full.names = TRUE)

# The location to print out the background buffers (.shp) (will be created if it doesn't exists)
buff_output <- "C:/Users/bshipley6/Desktop/EXAMPLE/TestRun/buffers"

# Generates buffers for each species.
# ncores should be set no higher than the number of cores the computer has minus 1
BackgroundBuffers(occlist = occlist,
                  envdata = TSEnv$training,
                  buff_output,
                  ncores = 2)

# 6. Background Points------------------------------

# This function creates the randomly generated background points and environmentally subsamples them
# (if required).

# Set the parameters for the background point generation (how many points, and how spatially-constrained)
nbg <- 1000
spatial_weights <- 0.5

# Make a list of the buffer files (generated in previous function)
bufflist <- list.files(buff_output, pattern = ".shp$", full.names = TRUE)

# Define the location where the background points will be printed out to (.shp) (will be created if it doesn't exist)
bg_output <- "C:/Users/bshipley6/Desktop/EXAMPLE/TestRun/backgrounds"
BackgroundPoints(spplist = spplist,
                 envdata = TSEnv$training,
                 output = bg_output,
                 nbg = nbg,
                 spatial_weights = spatial_weights,
                 buffers = bufflist,
                 method = "Varela",
                 ncores = 2)

# 7. Varying Environmental Sets by Species----------

# SDMs are generally better when a species-specific set of environmental variables are used (with biological relevance)
# This function can change the set of environmental variables used for each species to tailor the analysis to the
# habitat suitability of each species.

# Define a list of the environmental variables to keep for each species
envvar <- rep("Bio1,Bio12,Bio14,Bio6,Bio9", length = length(occlist))

# Define a list of the background point files
bglist <- list.files(bg_output, pattern = ".csv", full.names = TRUE)

# Like before, in this example, megaSDM overwrites the occurrence and background points,
# but they could be placed in a different folder if requested.
VariableEnv(occlist = occlist,
           bglist = bglist,
           env_vars = envvar,
           occ_output = occ_output,
           bg_output = bg_output)

# 8. MaxEnt Modelling--------------------------------

# After the occurrence and background point management,
# megaSDM will model the habitat suitability using the MaxEnt algorithm.
bglist <- list.files(bg_output, pattern = ".csv", full.names = TRUE)
# This is where the results of the MaxEnt model rusn will be printed out to (as .lambdas files)
# unlike the rest of the output folders, THIS ONE MUST BE CREATED BEFOREHAND
# It needs to have the maxent.jar file in it, so copy and paste the maxent.jar file (downlaod with this example)
# from its current location to the location given by "model_output"
model_output <- "C:/Users/bshipley6/Desktop/EXAMPLE/TestRun/models"

# This function generates the actual MaxEnt models for each species. The MaxEnt parameters
# (e.g., regualarization, crossvalidation meethod) can be manipluted within this function:
# see the documentation for a compelte list of settings that can be changed.

# "nrep" is to 4, meaning that the MaxEnt algorithm will run 4 times with a different subset of occurrence points
# for a better representation of the habitat suitability.
MaxEntModel(occlist = occlist,
            bglist = bglist,
            model_output = model_output,
            ncores = 2,
            nrep = 4,
            alloutputs = FALSE)


# 9. MaxEnt Projection------------------------------

# After generating the model parameters, the next function takes the outputted .lambdas files
# and projects them onto the current and future environmental conditions in the study region.

# First, create a list of the time periods and climate scenarios used in the analysis
# (starting with the year the model is trained on)
time_periods <- c(2010,2050,2070)
scenarios <- c("RCP4.5", "RCP8.5")

# Define the directory where the current study area rasters are located (written out in Step 1)
study_dir <- "C:/Users/bshipley6/Desktop/EXAMPLE/TestRun/studyarea"

# Define the directories where the future study area rasters are location (written out in Step 2)
# This argument should be a list of directories, separated into the different climate scenarios and years:
# list(c(Scenario1Year1, Scenario1Year2), c(Scenario2Year1, Scenario2Year2))

predictdir <- list(c("C:/Users/bshipley6/Desktop/EXAMPLE/TestRun/RCP4.5/2050",
                     "C:/Users/bshipley6/Desktop/EXAMPLE/TestRun/RCP4.5/2070"),
                   c("C:/Users/bshipley6/Desktop/EXAMPLE/TestRun/RCP8.5/2050",
                     "C:/Users/bshipley6/Desktop/EXAMPLE/TestRun/RCP8.5/2070"))

# Where the results will be printed out
result_dir <- "C:/Users/bshipley6/Desktop/EXAMPLE/TestRun/Results"

# This function takes the generated model parameters, projects them onto current and future environments,
# removes poorer models from the analyses, creates ensemble models from the multiple replicates conducted
# in the previous step, and generates binary presence/absence distribution maps.

# Other options are also available (check the documentation page)
MaxEntProj(input = model_output,
           time_periods = time_periods,
           scenarios = scenarios,
           study_dir = study_dir,
           predict_dirs = predictdir,
           output = result_dir,
           aucval = 0.7,
           ncores = 2)


# 10. TimeMaps---------------------------------------

# megaSDM can use the binary maps at each time period to create "time maps":  These
# maps detail the step-wise expansions and contractions of the species distribution
# through those time-steps, allowing for the visualization of both unidirectional
# range shifts (e.g., a range expansion across all time-steps) and more complex
# dynamics (e.g., a range expansion from 2010-2050 followed by a range contraction from 2050-2070).

# The time maps will be written out to the directory supplied in "result_dir"
createTimeMaps(result_dir = result_dir,
               time_periods = time_periods,
               scenarios = scenarios,
               dispersal = FALSE,
               ncores = 2)

# 11. AdditionalStats-------------------------------

# This function generates statistics on species range sizes and changes through the multiple time steps
# and different scenarios. It also craetes several graphs (written out to the directory provided by "result_dir")
# to visualize these changes.

additionalStats(result_dir = result_dir,
                time_periods = time_periods,
                scenarios = scenarios,
                dispersal = FALSE,
                ncores = 2)

# 12. dispersalRate---------------------------------

# Incorporating dispersal limitations into SDMs can help to mroe accurately predict where a species could live
# in the future. This function uses data on the average yearly dispersal rate of a species (given in km/year) to
# constrain the continuous habitat suitability maps and the presence/absence maps by the dispersal ability of each
# species.

# Add in dispersal data (normally you would read a .csv file with two columns, but I just added the data in by hand here)
dispersaldata <- data.frame(Species = spplist, Rate = c(48.92,
                                                        1.21,
                                                        6.07,
                                                        3.37,
                                                        0.96,
                                                        3.27))

dispersalRate(result_dir = result_dir,
              dispersaldata = dispersaldata,
              time_periods = time_periods,
              scenarios = scenarios,
              ncores = 2)

# 12. TimeMaps/AdditionalStats (with dispersal rate)-----------------------

# Repeat the time map and additional stats steps (Steps 10+11) for the dispersal constrained data
# dispersal = TRUE this time
# The additional stats function will compare the species ranges between the dispersal constrained adn the
# regular data
createTimeMaps(result_dir, time_periods, scenarios, dispersal = TRUE, dispersaldata = dispersaldata, ncores = 2)
additionalStats(result_dir, time_periods, scenarios, dispersal = TRUE, dispersaldata = dispersaldata, ncores = 2)

# 13. Species Richness Maps--------------------------

# Finally, this function creates species richness maps for the set of species examined across all
# time periods and scenarios. In addition, it compares the species richness of the dispersal-constrained
# analysis to the regular analysis.
createRichnessMaps(result_dir, time_periods, scenarios, dispersal = TRUE, taxonlist = FALSE)

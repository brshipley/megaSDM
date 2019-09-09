#Initializations--------------------------------
# Set location and name of config txt file to be used
configTXTloc <- paste(getwd(), "/config.txt", sep="")
print(paste("Current config file location: ", configTXTloc, sep=""))

# Pull data from config.txt file
config <- read.csv(configTXTloc, sep = "=", comment.char = "#", header = FALSE, stringsAsFactors = FALSE, strip.white = TRUE) 
row.names(config) <- config[, 1]
config <- subset(config, select = -V1)
names(config) <- "Parameters"

#Functions--------------------------------------
# Create subfolder and file locations based on main folder
createLocation <- function(mainFolder, dirloc) {
  return(paste(mainFolder, config[dirloc, ], sep = ""))
}

# Check if this experiment location exists; create if applicable
checkDirExists <- function(folderVariable) {
  if(!dir.exists(folderVariable)) {
    dir.create(folderVariable)
    print(paste("Directory Created: ", folderVariable, sep = ""))
  } 
}

# Check if file exists. Output to user to add file if applicable
# userInput = 1 if file has to be inputted by user; 0 if file is created by program
checkFileExists <- function(fileVariable, userInput) {
  # If user needs to input file but it isn't there
  if((userInput == 1) & (!file.exists(fileVariable))) {
    stop(paste("File doesn't exist: ", fileVariable, ". Please create the file and continue!"))
  }
}

# Check if file already exists, and give warning that it will be overwritten
fileExistsWarning <- function (fileVariable) {
  if (file.exists(fileVariable)) {
    message(paste("File already exists and may be overwritten: ", fileVariable, ". Please move or delete before you run!"))
  }
}

# Check if folder contains data, and give warning that it will be overwritten; maxStartFiles = max allowable folders/files to start
folderHasDataWarning <- function (folderVariable, maxStartFiles) {
  if (length(dir(folderVariable)) > maxStartFiles) {
    message(paste0("Folder already contains data and may be overwritten: ", folderVariable, ". Please move or delete files in these directories before you run!"))
  }
}

#Location Creation--------------------------
# Directories defined in config.txt created and added to df
TrialDirectory <- config["TrialDirectory", ]
spplist <- createLocation(TrialDirectory, "spplist")
df <- data.frame(spplist, stringsAsFactors = FALSE)
test <- createLocation(TrialDirectory, "test")
df <- data.frame(df, test, stringsAsFactors = FALSE)
samples <- createLocation(test, "samples")
df <- data.frame(df, samples, stringsAsFactors = FALSE)
occurrences <- createLocation(TrialDirectory, "occurrences")
df <- data.frame(df, occurrences, stringsAsFactors = FALSE)
buff_dir <- createLocation(TrialDirectory, "buff_dir")
df <- data.frame(df, buff_dir, stringsAsFactors = FALSE)
sppcountsloc <- createLocation(TrialDirectory, "sppcountsloc")
df <- data.frame(df, sppcountsloc, stringsAsFactors = FALSE)
result_dir <- createLocation(TrialDirectory, "result_dir")
df <- data.frame(df, result_dir, stringsAsFactors = FALSE)
counts <- createLocation(TrialDirectory, "counts")
df <- data.frame(df, counts, stringsAsFactors = FALSE)

#Run---------------------------------------
#Steps ("Y" or "N") added to df
ClipEnvDataStep <- config["ClipEnvDataStep", ]
df <- data.frame(df, ClipEnvDataStep, stringsAsFactors = FALSE)
CoordinateProjectionStep <- config["CoordinateProjectionStep", ]
df <- data.frame(df, CoordinateProjectionStep, stringsAsFactors = FALSE)
gbifstep <- config["gbifstep", ]
df <- data.frame(df, gbifstep, stringsAsFactors = FALSE)
OccurEnvFeaturesStep <- config["OccurEnvFeaturesStep", ]
df <- data.frame(df, OccurEnvFeaturesStep, stringsAsFactors = FALSE)
subsampleVarelaStep <- config["subsampleVarelaStep", ]
df <- data.frame(df, subsampleVarelaStep, stringsAsFactors = FALSE)
speciesBufferStep <- config["speciesBufferStep", ]
df <- data.frame(df, speciesBufferStep, stringsAsFactors = FALSE)
backgroundPointsStep <- config["backgroundPointsStep", ]
df <- data.frame(df, backgroundPointsStep, stringsAsFactors = FALSE)
UrbanAnalysis <- config["UrbanAnalysis", ]
df <- data.frame(df, UrbanAnalysis, stringsAsFactors = FALSE)
ProtectedAnalysis <- config["ProtectedAnalysis", ]
df <- data.frame(df, ProtectedAnalysis, stringsAsFactors = FALSE)
dispersalStep <- config["dispersalStep", ]
df <- data.frame(df, dispersalStep, stringsAsFactors = FALSE)
RichnessStep <- config["RichnessStep", ]
df <- data.frame(df, RichnessStep, stringsAsFactors = FALSE)
dispersalRan <- "N"
df <- data.frame(df, dispersalRan, stringsAsFactors = FALSE)

#Other preferences added to df
TrainingAreaClip <- config["TrainingAreaClip", ]
df <- data.frame(df, TrainingAreaClip, stringsAsFactors = FALSE)
TrainClipLatitude <- config["TrainClipLatitude", ]
df <- data.frame(df, TrainClipLatitude, stringsAsFactors = FALSE)
TrainClipLongitude <- config["TrainClipLongitude", ]
df <- data.frame(df, TrainClipLongitude, stringsAsFactors = FALSE)
decimalLatitude <- config["decimalLatitude", ]
df <- data.frame(df, decimalLatitude, stringsAsFactors = FALSE)
decimalLongitude <- config["decimalLongitude", ]
df <- data.frame(df, decimalLongitude, stringsAsFactors = FALSE)
minlat <- as.numeric(config["minlat_sa", ])
df <- data.frame(df, minlat, stringsAsFactors = FALSE)
maxlat <- as.numeric(config["maxlat_sa", ])
df <- data.frame(df, maxlat, stringsAsFactors = FALSE)
minlong <- as.numeric(config["minlong_sa", ])
df <- data.frame(df, minlong, stringsAsFactors = FALSE)
maxlong <- as.numeric(config["maxlong_sa", ])
df <- data.frame(df, maxlong, stringsAsFactors = FALSE)
rastertype <- config["rastertype", ]
df <- data.frame(df, rastertype, stringsAsFactors = FALSE)
ncores <- as.numeric(config["ncores", ])
df <- data.frame(df, ncores, stringsAsFactors = FALSE)
nsubsamp <- as.numeric(config["nsubsamp", ])
df <- data.frame(df, nsubsamp, stringsAsFactors = FALSE)
nclimatebins <- as.numeric(config["nclimatebins", ])
df <- data.frame(df, nclimatebins, stringsAsFactors = FALSE)
nPCAxes <- as.numeric(config["nPCAxes", ])
df <- data.frame(df, nPCAxes, stringsAsFactors = FALSE)
nrep <- as.numeric(config["nrep", ])
df <- data.frame(df, nrep, stringsAsFactors = FALSE)
reptype <- config["reptype", ]
df <- data.frame(df, reptype, stringsAsFactors = FALSE)
aucval <- as.numeric(config["aucval", ])
df <- data.frame(df, aucval, stringsAsFactors = FALSE)
hinge <- as.character(config["hinge", ])
df <- data.frame(df, hinge, stringsAsFactors = FALSE)
threshold <- as.character(paste0(config["threshold", ], ".logistic.threshold"))
df <- data.frame(df, threshold, stringsAsFactors = FALSE)
AllOutputs <- as.character(config["AllOutputs", ])
df <- data.frame(df, AllOutputs, stringsAsFactors = FALSE)
defaultCRS <- config["defaultCRS", ]
df <- data.frame(df, defaultCRS, stringsAsFactors = FALSE)
desiredCRS <- config["desiredCRS", ]
df <- data.frame(df, desiredCRS, stringsAsFactors = FALSE)
randomseed <- as.numeric(config["randomseed", ])
df <- data.frame(df, randomseed, stringsAsFactors = FALSE)
set.seed(randomseed)
currentyear <- as.numeric(config["currentyear", ])
df <- data.frame(df, currentyear, stringsAsFactors = FALSE)
resolution <- as.numeric(config["resolution", ])
df <- data.frame(df, resolution, stringsAsFactors = FALSE)
Categorical <- as.character(config["Categorical", ])
df <- data.frame(df, Categorical, stringsAsFactors = FALSE)
nbg <- as.character(config["nbg", ])
df <- data.frame(df, nbg, stringsAsFactors = FALSE)

#Checking the number of years
numYear <- length(grep("^Year", row.names(config)))
df <- data.frame(df, numYear, stringsAsFactors = FALSE)
for (i in 1:numYear) {
  yearName <- paste0("Year", i)
  Year <- as.numeric(config[yearName, 1])
  if (is.na(Year)) {
    stop("NumYear field in configuration file doesn't match the number of years! Please revise your configuration file!")
  } else {
    df <- data.frame(df, Year, stringsAsFactors = FALSE)
  }
}

years <-  rep(NA,length = numYear)
years <-  as.numeric(config[grep("^Year", row.names(config)), ])

#Checking the number of climate scenarios
numScenario <- length(grep("^Scenario", row.names(config)))
df <- data.frame(df, numScenario, stringsAsFactors = FALSE)
if (numScenario > 0) {
  for (i in 1:numScenario) {
    ScenName <- paste0("Scenario", i)
    Scenario <- config[ScenName, ]
    if (is.na(Scenario)) {
      stop("NumScenario field in configuration file doesn't match the number of scenarios! Please revise your configuration file!")
    } else {
      df <- data.frame(df, Scenario, stringsAsFactors = FALSE)
    }
  }
}


# All data used for the species distribution modelling, etc.
DataDirectory <- config["DataDirectory", ]
scripts <- createLocation(DataDirectory, "scripts")
df <- data.frame(df, scripts, stringsAsFactors = FALSE)
proj_trainingarea <- createLocation(DataDirectory, "proj_trainingarea")
df <- data.frame(df, proj_trainingarea, stringsAsFactors = FALSE)
trainingarea <- createLocation(DataDirectory, "trainingarea")
df <- data.frame(df, trainingarea, stringsAsFactors = FALSE)
proj_studyarea<- createLocation(DataDirectory, "proj_studyarea")
df <- data.frame(df, proj_studyarea, stringsAsFactors = FALSE)
studyarea<- createLocation(DataDirectory, "studyarea")
df <- data.frame(df, studyarea, stringsAsFactors = FALSE)
if (numScenario > 0) {
  proj_predictenv <- createLocation(DataDirectory, "proj_predictenv")
  df <- data.frame(df, proj_predictenv, stringsAsFactors = FALSE)
  predictenv <- createLocation(DataDirectory, "predictenv")
  df <- data.frame(df, predictenv, stringsAsFactors = FALSE)
}
if (ProtectedAnalysis == "Y") {
  protected_dir <- createLocation(DataDirectory, "protected_dir")
  df <- data.frame(df, protected_dir, stringsAsFactors = FALSE)
  proj_protected_dir <- createLocation(DataDirectory, "proj_protected_dir")
  df <- data.frame(df, proj_protected_dir, stringsAsFactors = FALSE)
}
if (UrbanAnalysis == "Y") {
  urbanized_dir <- createLocation(DataDirectory, "urbanized_dir")
  df <- data.frame(df, urbanized_dir, stringsAsFactors = FALSE)
  proj_urbanized_dir <- createLocation(DataDirectory, "proj_urbanized_dir")
  df <- data.frame(df, proj_urbanized_dir, stringsAsFactors = FALSE)
}
if (dispersalStep=="Y") {
  dispersalRate_dir <- createLocation(DataDirectory, "dispersalRate_dir")
  df <- data.frame(df, dispersalRate_dir, stringsAsFactors = FALSE)
}

# Check if folders exist
checkDirExists(TrialDirectory)
checkDirExists(test)
checkDirExists(sppcountsloc)
checkDirExists(occurrences)
checkDirExists(buff_dir)
checkDirExists(result_dir)

# Check if files already exist (warn that they will be overwritten)
if (gbifstep == "Y") {
  fileExistsWarning(counts)
  # Check if folder has data (warn that they will be overwritten)
  folderHasDataWarning(sppcountsloc, 0)
  folderHasDataWarning(test, 0) ## Test folder should start empty
}

return(df)
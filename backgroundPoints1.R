####backgroundPoints1.R####
##generates full-training area backgroundpoints (not species-specific)

#Initializations--------------------------------------
#Loads the necessary packages
library(rgeos)
library(raster)
library(dplyr)

#Loads the necessary variables from "df"
test <- df[, "test"]
proj <- df[, "proj_trainingarea"]
result_dir <- df[, "result_dir"]
buff_dir <- df[, "buff_dir"]
nsubsamp <- as.numeric(df[, "nsubsamp"])
samples <- "samples"
nPCAxes <- df[, "nPCAxes"]
rastertype<- df[, "rastertype"]
desiredCRS <- df[, "desiredCRS"]
randomseed <- df[, "randomseed"]
Categorical <- df[, "Categorical"]
spatial_weights <- df[, "spatial_weights"]

ncores <- df[, "ncores"]
outfile <- "statsout.txt"
speciesBufferStep <- df[, "speciesBufferStep"]
spp_total <- spp_total

nbg <- as.numeric(df[, "nbg"])
if ((speciesBufferStep == "Y") | (length(list.files(buff_dir, pattern = paste0(rastertype, "$|.shp"))) == length(spp_total))) {
  nbgTrain <- nbg * (1 - spatial_weights)
} else {
  nbgTrain <- nbg
}

#To produce the same background csv files every time the code is run with the same inputs
if (df[, "randomseed"]) {
  set.seed(randomseed)
}
setwd(test)

#creates a new folder "backgrounds" in "test"
dir.create("backgrounds")

#disables scientific notation
options(scipen = 999)

#creates a stack of all environmental rasters as found in "proj"
#the stack is called "train" for "Training Area"
train <- list.files(path = proj, full.names = TRUE, pattern = paste0("\\.bil$"))
train <- stack(train)
crs(train) <- desiredCRS

#If nbg is more than the number of cells in the training raster, nbg should be equal to the number of cells
if (nbgTrain > ncell(train[[1]])) {
  nbgTrain <- ncell(train[[1]])
  message(paste0("Warning: The number of background points is more than the number of raster cells (", ncell(train[[1]])), ")")
  message(paste0("    Changing the number of background points to ", ncell(train[[1]])))
}


#Functions-----------------------------------------
#VarelaSample does environmental subsampling of the generated background points (10*number of background points desired)
VarelaSample <- function (occurrences, ClimOccur, no_bins) {
  nsamples <- c()
  for (j in 1:length(no_bins)) {
    out_ptz <- EnvOccur[, 1:2]
    for (i in 3:length(names(EnvOccur))) {
      #make a data frame that is this variable with no NA values
      k <- EnvOccur[!is.na(EnvOccur[, i]), i]
      #calculate the observed range of this variable
      rg <- range(k)
      #figure out the resolution from the number of bins
      res <- (rg[2] - rg[1]) / no_bins[j]
      #rescale the axis by the range and bin size, so the value is just a 
      #number from 1 to no_bins for its bin membership
      d <- (EnvOccur[, i] - rg[1]) / res
      #d is now a vector of values ranging from 0 to no_bins
      f <- ceiling(d)
      #f is a vector of bin membership
      f[f == 0] <- 1 #move the zeros into the 1 bin
      #correct the name of the vector, so it will carry over to the output
      names(f) <- names(EnvOccur)[i]
      #add the bin membership vector to the output df for this section
      out_ptz <- cbind(out_ptz, f)
      #get the names correct
      names(out_ptz)[length(names(out_ptz))] <- names(EnvOccur)[i]
    }
    
    #subsample the bin membership df to come up with the filled bins
    sub_ptz <- distinct(out_ptz[, -1:-2])
    #count the number of filled bins
    no_grps <- nrow(sub_ptz)
    #add a column with the group membership number; this number is arbitrary
    sub_ptz$grp <- c(1:no_grps)
    
    #join the out_ptz with the subsample to capture the group membership info
    #note: join() will automatically match the variable names from these two dfs
    out_ptz <- left_join(out_ptz, sub_ptz)
    #out_ptz now has a group membership  for each input point
    
    #select a random point for each group -- this is an additional improvement on the 
    #Varela et al. function, because it does not pick the same points each time.
    
    #make a landing spot for the data
    final_out <- data.frame(x = numeric(), y = numeric())
    
    #cycle through each group
    for (i in 1:no_grps) {
      #subset to the members of the ith group, keep only the Latitude and Longitude
      grp_mbrs <- out_ptz[out_ptz$grp == i, c(1, 2)]
      
      #pick one of these group members to output
      grp_out <- grp_mbrs[sample(1:nrow(grp_mbrs), 1), ]
      #bind this sampled row to the output df
      final_out <- rbind(final_out, grp_out)
    }
    
    #return the subsampled points as a df of Latitude and Longitude values
    final_out <- data.frame(x = final_out[, 1],y = final_out[, 2])
    final_out <- merge(final_out, ClimOccur, by = c("x","y"), all.x = TRUE)
    nsamples <- c(nsamples, nrow(final_out))
  }
  if (length(no_bins) == 1) {
    return(final_out)
  } else {
    return(data.frame(NumberofSamples = nsamples,NumberOfClimateBins = no_bins))
  }
}

#Run-----------------------------------------------
#Generates empty vectors for background results 
BuffClimateBins <- c()
BuffNumber <- c()
TrainClimateBins <- c()
TrainNumber <- c()
TotalNumber <- c()

if (nbgTrain > 0) {
  for (i in 1:nsubsamp) {
    #Extract many more climate points than we need (will get subsampled down to correct number later)
    #if nbg*10 is more than the number of cells in the raster, then sampleRandom won't work
    if ((nbgTrain * 10) > ncell(train[[1]])) {
      RastTrainBG <- ncell(train[[1]])
    } else {
      RastTrainBG <- nbgTrain * 10
    }
    ClimOccur <- sampleRandom(train, RastTrainBG, na.rm = TRUE, xy = TRUE)
  
    #Run PCA analysis
    if (is.na(Categorical)) {
      PCAEnv <- prcomp(ClimOccur[, 3:ncol(ClimOccur)], scale = TRUE)
      PCAImp <- summary(PCAEnv)$importance
      
      #Determine the number of PC axes to use for subsampling
      if (!is.na(nPCAxes)) {
        NumberAxes <- nPCAxes
      } else {
        NumberAxes <- max(2, min(which(PCAImp[3,] > 0.95)))
      }
      
      if (NumberAxes > nlayers(train)) {
        NumberAxes <- nlayers(train)
        message(paste0("'NPCAxes' is more than the number of environmental variables (", NumberAxes, "): Setting 'NPCAxes' to ", NumberAxes, "."))
      }
      
      #Add PCA values to the unsubsampled data frame
      EnvOccur <- data.frame(cbind(ClimOccur[, 1:2], PCAEnv$x[, 1:NumberAxes]))
    } else {
      EnvOccur <- as.data.frame(ClimOccur)
    }
    
    #Determine the number of climate bins that are needed to make nbgTrain points
    #For the first run, start at 1. Every other time, start at nbinscount - 10
    if (i == 1) {
      nbinscount <- 1
    } else {
      nbinscount <- max(1, nbinscount - 10)
    }
    
    #Steadily increase the number of bins until the desired nbg is reached
    nbgTest <- 0
    bgdf_train <- data.frame(c())
    while (nbgTest < nbgTrain && nbinscount < 99) {
      bgdf_train <- VarelaSample(EnvOccur[, 1:2], ClimOccur, c(nbinscount:(nbinscount+1)))
      nbinscount <- nbinscount + 1
      nbgTest <- bgdf_train[nrow(bgdf_train), 1]
    }
    
    #Output the number of background points and climate bins
    NBinsFinal <- bgdf_train[which(abs(bgdf_train[, 1] - nbgTrain) == min(abs(bgdf_train[, 1] - nbgTrain)))[1], 2]
    print(paste0("Number of Climate Bins Used For Training Area Subsampling: ", NBinsFinal))
    
    bgdf_train <- VarelaSample(RandomTrain[, 1:2], ClimOccur, NBinsFinal)
    print(paste0("Number of Background Points In Training Area: ", nrow(bgdf_train)))
    
    #Create data frame of background points
    bgtraindataframe <- data.frame(bgdf_train, stringsAsFactors = FALSE)
    
    #Statistics about background points/number of bins
    TrainClimateBins <- c(TrainClimateBins, NBinsFinal)
    TrainNumber <- c(TrainNumber, nrow(bgtraindataframe))
    
    #Write background CSV Files
    write.csv(bgtraindataframe, row.names = FALSE, file = paste0("backgrounds/Train_Background_", i, ".csv"))
  }
} else {
  for (i in 1:nsubsamp) {
    bgtraindataframe <- sampleRandom(train, 1, na.rm = TRUE, xy = TRUE)
    bgtraindataframe <- bgtraindataframe[-1, ]
    #Write background CSV Files
    write.csv(bgtraindataframe, row.names = FALSE, file = paste0("backgrounds/Train_Background_", i, ".csv"))
  }
  TrainClimateBins <- rep(0, nsubsamp)
  TrainNumber <- rep(0, nsubsamp)
}
#Write Statistics CSV into Results Folder
bgstats <- data.frame(Subsample = c(1:nsubsamp), TrainClimateBins, TrainNumber)

for (i in nrow(bgstats)) {
  if (bgstats$TrainClimateBins[i] == 99) {
    message("Warning: the number of background points wanted is too high! Only 99 environmental bins used")
  }
}
write.csv(bgstats, file = paste0(result_dir, "/BackgroundPoints_stats.csv"), row.names = FALSE)
rm(train, RastTrainBG)
gc()
####backgroundPoints2.R####
##Generates buffer area backgroundpoints and/or species-specific background points

#Initializations----------------------------------------------
#Loads the necessary packages
library(rgeos)
library(raster)
library(dplyr)

#Loads the necessary variables from "df"
test <- df[, "test"]
proj <- df[, "proj_trainingarea"]
buff_dir <- df[, "buff_dir"]
result_dir <- df[, "result_dir"]
nsubsamp <- as.numeric(df[, "nsubsamp"])
samples <- "samples"
nbins <- df[, "nclimatebins"]
nPCAxes <- df[, "nPCAxes"]
Categorical <- df[, "Categorical"]
spatial_weights <- df[, "spatial_weights"]
rastertype<- df[, "rastertype"]
desiredCRS <- df[, "desiredCRS"]
randomseed <- df[, "randomseed"]
ncores <- df[, "ncores"]
nbg <- df[, "nbg"]
outfile <- "statsout.txt"
speciesBufferStep <- df[, "speciesBufferStep"]

if (df[, "randomseed"]) {
  set.seed(randomseed)
}

#Generates a global bounding box to remove problems with the anti-meridian
GlobalBBox <- SpatialPoints(matrix(c(-179.9, 179.9, 89.9, -89.9), nrow = 2), proj4string = crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs"))
setwd(test)

#creates a new folder "backgrounds" in "test"
if (!dir.exists(paste0(test, "/backgrounds"))) {
  dir.create("backgrounds")
}

#Creates a list of the active species including the "samples" path
ListSpp <- list.dirs(path = samples, recursive = FALSE)

#disables scientific notation
options(scipen = 999)

#creates a stack of all environmental rasters as found in "proj"
#the stack is called "train" for "Training Area"
train <- list.files(path = proj, full.names = TRUE, pattern = paste0("\\.bil$"))
train <- stack(train)
crs(train) <- desiredCRS

#Functions----------------------------------
#This function environmentally subsamples the background points generated
VarelaSample <- function (OccurData, ClimOccur, no_bins) {
  
  #make a landing spot for the bin membership vectors
  #cycle through all of the environmental variables (columns 3 to end)
  nsamples <- c()
  for (j in 1:length(no_bins)) {
    out_ptz <- OccurData[, 1:2]
    for(i in 3:ncol(OccurData)) {
      #make a data frame that is this variable with no NA values
      k <- OccurData[!is.na(OccurData[, i]), i]
      #calculate the observed range of this variable
      rg <- range(k)
      #figure out the resolution from the number of bins
      res <- (rg[2] - rg[1]) / no_bins[j]
      #rescale the axis by the range and bin size, so the value is just a 
      #number from 1 to no_bins for its bin membership
      d <- (OccurData[, i] - rg[1]) / res
      #d is now a vector of values ranging from 0 to no_bins
      f <- ceiling(d)
      #f is a vector of bin membership
      f[f == 0] <- 1 #move the zeros into the 1 bin
      #correct the name of the vector, so it will carry over to the output
      names(f) <- names(OccurData)[i]
      #add the bin membership vector to the output df for this section
      out_ptz <- cbind(out_ptz, f)
      #get the names correct
      names(out_ptz)[length(names(out_ptz))] <- names(OccurData)[i]
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
    for(i in 1:no_grps) {
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
    return(data.frame(NumberofSamples = nsamples, NumberOfClimateBins = no_bins))
  }
}

bgpoints <- function(CurSpp) {
  
  #Generates empty vectors for background results
  BuffNumber <- c()
  TrainNumber <- c()
  TotalNumber <- c()
  
  #creates an active species, "spp" and then cuts it to only the spp name
  spp <- CurSpp
  spp.name <- substr(spp, 9, (nchar(spp) - 4))
  
  #If species-dependent numbers of background points are necessary
  if (length(grep("x", tolower(nbg))) > 0) {
    #Calculate number of occurrences
    OccurrenceFileName <- list.files(path = paste0(test, "/", samples, "/", spp.name), pattern = paste0("\\.csv$"))[1]
    OccurrenceFile <- read.csv(paste0(test, "/", samples, "/", spp.name, "/", OccurrenceFileName))
    OccurrenceNumber <- nrow(OccurrenceFile)
    
    #Multiply by the desired number
    Multiplier <- as.numeric(substr(nbg, 1, nchar(nbg) - 1))
    nbg <- Multiplier * OccurrenceNumber
    if (nbg > ncell(train[[1]])) {
      nbg <- ncell(train[[1]])
      message(paste0("Warning: The number of background points is more than the number of raster cells (", ncell(train[[1]])), ")")
      message(paste0("    Changing the number of background points to ", ncell(train[[1]])))
    }
  } else {
    nbg <- as.numeric(nbg)
    if (nbg > ncell(train[[1]])) {
      nbg <- ncell(train[[1]])
      message(paste0("Warning: The number of background points is more than the number of raster cells (", ncell(train[[1]])), ")")
      message(paste0("    Changing the number of background points to ", ncell(train[[1]])))
    }
  }
  
  #For each desired occurrence subsample
  for (g in 1:nsubsamp) {
    #Step 1: sampling the points inside the buffer (if exists)----------
    setwd(buff_dir)
    
    if (length(list.files(buff_dir)) > 0) {
      #Calculate the number of buffer points necessary
      nbgBuffer <- nbg * spatial_weights
      
      if(round(nbgBuffer) == 0) {
        nbgBuffer <- 1
      }
      
      #if buffers are rasters, use CRS. Otherwise use SpatialPolygons
      if (file.exists(paste0(buff_dir, "/", spp.name, rastertype))) {
        BuffFile <- raster(paste0(spp.name, rastertype), native = TRUE)
        crs(BuffFile) <- desiredCRS
      } else {
        BuffFile <- shapefile(paste0(spp.name, ".shp"))
        if (file.exists(paste0(spp.name, ".prj"))) {
          maskTemp <- SpatialPolygons(BuffFile@polygons, proj4string = BuffFile@proj4string)
          GlobalBox2 <- spTransform(GlobalBBox, CRSobj = maskTemp@proj4string)
          maskTemp <- crop(maskTemp, extent(GlobalBox2))
          BuffFile <- spTransform(maskTemp, CRSobj = crs(desiredCRS))
        } else {
          BuffFile <- SpatialPolygons(BuffFile@polygons, proj4string = crs(desiredCRS))
          GlobalBox2 <- spTransform(GlobalBBox, CRSobj = crs(desiredCRS))
          BuffFile <- crop(BuffFile, extent(GlobalBox2))
        }
      }
      
      #If buffers are in raster form (else on line 260) 
      if (file.exists(paste0(spp.name, rastertype))) {
        #randomly sample %nbg*5% number of points from the buffer raster
        #if nbg*5 is more than the number of cells in the buffer raster, then simply sample all of the raster cells
        if ((nbg * 5) > cellStats(BuffFile, stat = sum)) {
          RastBuffBG <- cellStats(BuffFile, stat = sum)
        } else {
          RastBuffBG <- nbg * 5
        }
        
        #Randomly sample nbg * 5 points from the raster
        rasterValues <- rasterToPoints(BuffFile)
        randomCoord <- rasterValues[sample(nrow(rasterValues), size = RastBuffBG), ]
        randomCoordinates <- cbind(randomCoord[, 1], randomCoord[, 2])
        
        #generate SpatialPoints object from sampled raster values
        randomCoordinates <- SpatialPoints(randomCoordinates, proj4string = crs(desiredCRS))
        
        #Extract climate data
        ClimOccur <- extract(train, randomCoordinates)
        ClimOccur <- data.frame(x = randomCoordinates@coords[, 1], y = randomCoordinates@coords[, 2], ClimOccur)
        ClimOccur <- na.omit(ClimOccur)
        
        #Run PCA analysis
        if (is.na(Categorical) & round(nbg * spatial_weights) > 0) {
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
          EnvOccur <- cbind(ClimOccur[, 1:2], PCAEnv$x[, 1:NumberAxes])
        } else {
          EnvOccur <- ClimOccur
        }
        
        #Varela Subsamples the buffer background points with a range of climate bins to get the closest to nbg * spatial_weights
        nbinscount <- 1
        nbgTest <- 0
        bgdf_buffer <- data.frame(c())
        NUniqueClim <- nrow(unique(EnvOccur[, 3:ncol(EnvOccur)]))
        if (NUniqueClim < nbgBuffer) {
          nbgBuffer <- NUniqueClim
          message("Warning! The number of unique climates within the buffer is less than the number of background points wanted")
          message(paste0("Only ", NUniqueClim, " background points within the buffer can be generated"))
          nbinscount <- 98
        }
        
        #Steadily increase the number of bins until the desired nbg is reached
        while (nbgTest < nbgBuffer && nbinscount < max(99, nbins)) {
          bgdf_buffer <- VarelaSample(EnvOccur, ClimOccur, c(nbinscount:(nbinscount + 1)))
          nbinscount <- nbinscount + 1
          nbgTest <- bgdf_buffer[nrow(bgdf_buffer), 1]
        }
        
        #Output the number of background points and climate bins
        NBinsFinal <- bgdf_buffer[which(abs(bgdf_buffer[, 1] - nbgBuffer) == min(abs(bgdf_buffer[, 1] - nbgBuffer)))[1], 2]
        print(paste0("Number of Climate Bins Used For Buffer Subsampling: ", NBinsFinal))
        
        bgdf_buffer <- VarelaSample(EnvOccur, ClimOccur, NBinsFinal)
        print(paste0("Number of Background Points In Buffer: ", nrow(bgdf_buffer)))
      } else {
        #If the buffers are shapefiles
        
        #randomly samples nbg*5 number of points from the shapefile
        randomCoord <- spsample(BuffFile, (nbgBuffer * 5), "random", iter = 15)
        randomCoord2 <- spTransform(randomCoord, crs(train))
        randomCoordinates <- randomCoord2@coords
        
        #Extracts climate data
        ClimOccur <- extract(train, randomCoordinates)
        ClimOccur <- data.frame(x = randomCoord2@coords[, 1], y = randomCoord2@coords[, 2], ClimOccur)
        ClimOccur <- na.omit(ClimOccur)
        
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
          }
          
          #Add PCA values to the unsubsampled data frame
          EnvOccur <- cbind(ClimOccur[, 1:2], PCAEnv$x[, 1:NumberAxes])
        } else {
          EnvOccur <- ClimOccur
        }
        #Varela Subsamples the buffer background points with a range of climate bins to get the closest to nbg * spatial_weights
        nbinscount <- 1
        nbgTest <- 0
        bgdf_buffer <- data.frame(c())
        NUniqueClim <- nrow(unique(EnvOccur[, 3:ncol(EnvOccur)]))
        if (NUniqueClim < nbgBuffer) {
          nbgBuffer <- NUniqueClim
          message("Warning! The number of unique climates within the buffer is less than the number of background points wanted")
          message(paste0("Only ", NUniqueClim, " background points within the buffer can be generated"))
          nbinscount <- 98
        }
        
        #Steadily increase the number of bins until the desired nbg is reached
        while (nbgTest < nbgBuffer && nbinscount < max(99, nbins)) {
          bgdf_buffer <- VarelaSample(EnvOccur, ClimOccur, c(nbinscount:(nbinscount+1)))
          nbinscount <- nbinscount + 1
          nbgTest <- bgdf_buffer[nrow(bgdf_buffer), 1]
        }
        
        #Output the number of background points and climate bins
        NBinsFinal <- bgdf_buffer[which(abs(bgdf_buffer[, 1] - nbgBuffer) == min(abs(bgdf_buffer[, 1] - nbgBuffer)))[1], 2]
        print(paste0("Number of Climate Bins Used For Buffer Subsampling: ", NBinsFinal))
        
        bgdf_buffer <- VarelaSample(EnvOccur, ClimOccur, NBinsFinal)
        print(paste0("Number of Background Points In Buffer: ", nrow(bgdf_buffer)))
      }
      
      #adds species names & coordinates & makes random points+env data into a dataframe
      setwd(test)
      
      bgbuffdataframe <- data.frame(Species = rep(spp.name, nrow(bgdf_buffer)), bgdf_buffer)
      
      #removes rows containing NA
      bgbuffdataframe <- bgbuffdataframe[complete.cases(bgbuffdataframe), ]
      rowlength <- nrow(bgbuffdataframe)
      
      #Adds information to the stats output
      BuffNumber <- c(BuffNumber, rowlength)
      
      dir.create(paste0("backgrounds/", spp.name))
    } else {
      message ("No buffer file(s) found, skipping spatially constrained background point selection")
      bgbuffdataframe <- c()
      BuffNumber <- c(BuffNumber, "NA")
      dir.create(paste0(test, "/backgrounds/", spp.name))
    }
    
    #Step 2: Sampling from the entire training area------------ 
    #Checks to see if there are already background files created from "backgroundPoints1.R"
    BGTrainFiles <- list.files(path = paste0(test,"/backgrounds"), pattern = paste0("Train_Background_"))
    
    if (length(BGTrainFiles) != nsubsamp) {  
      #if nbg*5 is more than the number of cells in the buffer raster, then sampleRandom won't work
      if ((nbg * 5) > ncell(train[[1]])) {
        RastTrainBG <- ncell(train[[1]])
      } else {
        RastTrainBG <- nbg * 5
      }
      
      #Sample the required number of full-training area background points
      RandomTrain <- sampleRandom(train, RastTrainBG, na.rm=TRUE, xy=TRUE)
      RandomTrain <- as.data.frame(RandomTrain)
      
      #Run PC analysis
      if (is.na(Categorical)) {
        PCAEnv <- prcomp(RandomTrain[, 3:ncol(RandomTrain)], scale = TRUE)
        PCAImp <- summary(PCAEnv)$importance
        #Determine the number of PC axes to use for subsampling
        if (!is.na(nPCAxes)) {
          NumberAxes <- nPCAxes
        } else {
          NumberAxes <- max(2, min(which(PCAImp[3,] > 0.95)))
        }
        
        if (NumberAxes > nlayers(train)) {
          NumberAxes <- nlayers(train)
        }
        
        #Add PCA values to the unsubsampled data frame
        EnvOccur <- cbind(RandomTrain[, 1:2], PCAEnv$x[, 1:NumberAxes])
      } else {
        EnvOccur <- RandomTrain
      }
      #Gets the number of bins that will most closely approximate occurrences*10 (combined with sampled buffer points)
      if (length(list.files(buff_dir)) > 0) {
        nbinscount <- 1
        nbgTrain <- 0
        bgdf_train <- data.frame(c())
        
        #Steadily increase the number of bins until the desired nbg is reached
        while (nbgTrain < (nbg - nrow(bgbuffdataframe)) && nbinscount < max(99, nbins)) {
          bgdf_train <- VarelaSample(EnvOccur, RandomTrain, c(nbinscount:(nbinscount + 1)))
          nbinscount <- nbinscount + 1
          nbgTrain <- bgdf_train[nrow(bgdf_train), 1]
        }
        
        #Shorten by splitting up
        NBinsFinal <- bgdf_train[which(abs((bgdf_train$NumberofSamples + nrow(bgbuffdataframe)) - nbg) == min(abs(bgdf_train$NumberofSamples + nrow(bgbuffdataframe) - nbg)))[1], 2]
        if (length(NBinsFinal) > 1) {
          NBinsFinal <- NBinsFinal[1]
        }
        bgdf_train <- VarelaSample(EnvOccur, RandomTrain, c(NBinsFinal))
      } else {
        #if no buffer files exists, all points sampled from training area
        nbinscount <- 1
        nbgTrain <- 0
        bgdf_train <- data.frame(c())
        
        #Steadily increase the number of bins until the desired nbg is reached
        while (nbgTrain < nbg && nbinscount < max(99, nbins)) {
          bgdf_train <- VarelaSample(EnvOccur, RandomTrain, c(nbinscount:(nbinscount + 1)))
          nbinscount <- nbinscount + 1
          nbgTrain <- bgdf_train[nrow(bgdf_train), 1]
        }
        
        NBinsFinal <- bgdf_train[which(abs(bgdf_train$NumberofSamples - nbg) == min(abs(bgdf_train$NumberofSamples - nbg)))[1], 2]
        if (length(NBinsFinal) > 1) {
          NBinsFinal <- NBinsFinal[1]
        }
        bgdf_train <- VarelaSample(EnvOccur, RandomTrain, NBinsFinal)    
      }
      
      print(paste0("Number of Climate Bins Used For Training Area Subsampling: ", NBinsFinal))
      print(paste0("Number of Background Points In Training Area: ", nrow(bgdf_train)))
      
      #Creates the final background points data frame and append the statistics to the statistics vectors
      bgtraindataframe <- data.frame(Species = rep(spp.name,nrow(bgdf_train)), bgdf_train, stringsAsFactors = FALSE)
      
      if (length(list.files(buff_dir)) > 0 & (round(nbg* spatial_weights) != 0)) {
        Full_BGPoints <- rbind(bgbuffdataframe, bgtraindataframe)
      } else {
        Full_BGPoints <- bgtraindataframe
      }
      
      TrainNumber <- c(TrainNumber, nrow(bgtraindataframe))
      TotalNumber <- c(TotalNumber, nrow(Full_BGPoints))
      write.csv(Full_BGPoints, row.names = FALSE, file = paste0(test, "/backgrounds/", spp.name, "/", spp.name, "_background_", g, ".csv"))
    } else {
      
      print("Appending pre-generated training area background points...")
      bgdf_train <- read.csv(paste0(test, "/backgrounds/Train_Background_", g, ".csv"))
      bgtraindataframe <- data.frame(Species = rep(spp.name, nrow(bgdf_train)), bgdf_train)
      Full_BGPoints <- rbind(bgbuffdataframe, bgtraindataframe)
      TotalNumber <- c(TotalNumber, nrow(Full_BGPoints))
      write.csv(Full_BGPoints, row.names = FALSE, file = paste0(test, "/backgrounds/", spp.name, "/", spp.name, "_background_", g, ".csv"))
    }
  }
  #Creates (or appends) stats file for background points
  if (length(BGTrainFiles) == nsubsamp) {
    #If backgroundPoints1.R already generated a stats file, append new data
    BG_Stats <- data.frame(read.csv(paste0(result_dir, "/BackgroundPoints_stats.csv")))
    if (file.exists(paste0(result_dir, "/", spp.name, "/BackgroundPoints_stats.csv"))){
      Buffer_Stats <- data.frame(read.csv(paste0(result_dir, "/", spp.name, "/BackgroundPoints_stats.csv")))
      BG_Stats$BufferWidth <- Buffer_Stats$x
    }
    BG_Stats$BuffNumber <- BuffNumber
    BG_Stats$TotalNumber <- TotalNumber
    write.csv(BG_Stats, file = paste0(result_dir, "/", spp.name, "/BackgroundPoints_stats.csv"), row.names = FALSE)
  } else {
    #Creates statistics data frame
    BG_Stats <- data.frame(Subsample = c(1:nsubsamp), BuffNumber, TrainNumber, TotalNumber)
    if (file.exists(paste0(result_dir, "/", spp.name, "/BackgroundPoints_stats.csv"))){
      Buffer_Stats <- data.frame(read.csv(paste0(result_dir, "/", spp.name, "/BackgroundPoints_stats.csv")))
      BG_Stats$BufferWidth <- Buffer_Stats$x
    }
    write.csv(BG_Stats, file = paste0(result_dir, "/", spp.name, "/BackgroundPoints_stats.csv"), row.names = FALSE)
  }
  rm(BuffFile)
  gc()
}

#Run---------------------------
#Creates species list for parallelization
ListSpp <- c()
setwd(test)
speciesWorked <- spp_batch
for (i in 1:length(speciesWorked)) {
  ListSpp <- c(ListSpp, list.files(path = samples, full.names = TRUE, pattern = speciesWorked[i]))
}
ListSpp <- unique(ListSpp)
print("    Generating Background Points for:")
print(paste0("        ", spp_batch))

#Parallelization
clus <- makeCluster(ncores, outfile = outfile, setup_timeout = 0.5)
clusterExport(clus, varlist = c("test", "ncores", "proj", "buff_dir", "nsubsamp", "samples", "nbg",
                                "rastertype", "desiredCRS", "randomseed", "ListSpp", "speciesBufferStep",
                                "train", "bgpoints", "nPCAxes", "VarelaSample", "spatial_weights", "nbins", "result_dir", "GlobalBBox", "Categorical"))



clusterEvalQ(clus, library(rgeos))
clusterEvalQ(clus, library(raster))
clusterEvalQ(clus, library(dplyr))


print("    Creating background files in:")
print(paste0("        ", test, "/backgrounds/%Species_Name%"))

out <- parLapply(clus, ListSpp, function(x) bgpoints(x))
stopCluster(clus)
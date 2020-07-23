####UniqueSpeciesClimates.R####
##Limits background/occurrence point data to desired environmental variables
##UniqueSpeciesClimateStep = "Y"
##Use when modelling many species that do not have similar physiological requirements or constraints

#Initializations------------------------------
#Loads the necessary packages
library(parallel)

#Reads the species list (with the desired environmental variables for each species)
spplist <- df[, "spplist"]
spplist <- read.csv(spplist)

ncores <- as.numeric(df[, "ncores"])

test <- df[, "test"]
setwd(test)

#Functions------------------------------------------
#Removes uneeded environmental variables from the background and occurrence points
UniqueSpecClim <- function(CurSpp) {
  #Make lists of background and occurrence files for the current species
  BackgroundFiles <- list.files(path = paste0(getwd(), "/backgrounds/", CurSpp), full.names = TRUE)
  OccurrenceFiles <- list.files(path = paste0(getwd(), "/samples/", CurSpp), full.names = TRUE)
  
  #Determine which environmental variables were requested for that species
  ReqVarList <- as.vector(spplist[grep(paste0(gsub("_", " ", CurSpp), "$"), spplist[, 2]), 3])
  ReqVarList <- unlist(strsplit(ReqVarList, ",\\s*", perl = TRUE))
  ReqVar <- paste(paste0("^", ReqVarList, "$"), collapse = "|")
  
  #Subset the background files to only include the required environmental variables
  for (b in 1:length(BackgroundFiles)) {
    BGFile <- read.csv(BackgroundFiles[b])
    BGFileCoords <- BGFile[, c("Species", "x", "y")]
    ReqVarCol <- grep(tolower(ReqVar), tolower(names(BGFile)))
    if (length(ReqVarCol) != length(ReqVarList)) {
      stop(paste0("One or more climate variables given in 'spplist.csv' do not correspond to provided environmental rasters", "\n",  
                  "          Check the environmental variable column in 'spplist.csv'"))
    }
    BGFile_Final <- data.frame(BGFileCoords, BGFile[, ReqVarCol])
    write.csv(BGFile_Final, BackgroundFiles[b], row.names = FALSE)  
  }
  
  #Subset the occurrence files to only include the required environmental variables
  for (s in 1:length(OccurrenceFiles)) {
    OCFile <- read.csv(OccurrenceFiles[s])
    OCFileCoords <- OCFile[, c("Species", "x", "y")]
    ReqVarCol <- grep(tolower(ReqVar), tolower(names(OCFile)))
    if (length(ReqVarCol) != length(ReqVarList)) {
      stop("One or more climate variables given in 'spplist.csv' do not correspond to provided environmental rasters")
    }
    OCFile_Final <- data.frame(OCFileCoords, OCFile[, ReqVarCol])
    write.csv(OCFile_Final, OccurrenceFiles[s], row.names = FALSE)  
  }
}

#Run-------------------------------------
#Gets a list of the species analyzed on this go-around
ListSpp <- c()
speciesWorked <- spp_batch
for (i in 1:length(speciesWorked)) {
  ListSpp[i] <- paste0(substr(speciesWorked[i], 1, nchar(speciesWorked[i]) - 4))
}
ListSpp <- unique(ListSpp)
print("    Will evaluate species:")
print(paste0("        ", spp_batch))

#Parallelization
clus <- makeCluster(ncores, outfile = outfile, setup_timeout = 0.5)

clusterExport(clus, varlist = c("spplist", "test", "ListSpp", "UniqueSpecClim"))
out<-parLapply(clus, ListSpp, function(x) UniqueSpecClim(x))
stopCluster(clus)
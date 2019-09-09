#Initializations-------------------------------------
library(biomod2)
library(gtools)
library(parallel)
library(rgdal)
library(raster)
test <- df[, "test"]
ncores <- as.numeric(df[, "ncores"])
outfile <- "statsout.txt"
rastertype <- df[, "rastertype"]
desiredCRS <- df[, "desiredCRS"]
defaultCRS <- df[, "defaultCRS"]
dispersalStep <- df[, "dispersalStep"]
backgroundPointsStep <- df[, "backgroundPointsStep"]

setwd(df[, "proj_trainingarea"])
#load climate data:
train <- list.files(pattern = paste0('.bil$'), full.names = TRUE)
# sorting so they are in ascending order
train <- mixedsort(train)

train2 <- c()
for (i in 1:length(train)) {
  train2 <- c(train2, raster(train[i]))
}

train2 <- stack(train2)
crs(train2) <- desiredCRS

print("   Climate data loaded:")
print(paste0("      ", names(train2)))

#Functions---------------------------------------------------
run = function(CurSpp) {
  
  #read in the occurrence data frame (should be formatted with columns 'Species', 'Lat', 'Long')
  CurSpp2 <- (read.csv(CurSpp, stringsAsFactors = FALSE))
  
  #just in case other data came in... Also will ensure that data is in proper format
  
  tryCatch(CurSpp2 <- CurSpp2[, c("Species", "Longitude", "Latitude")],
           error = function(e) "Not the right column headings"
           )
  
  Species <- CurSpp2$Species
  
  #Convert data into a SpatialPoints format
  
  coordinates <- CurSpp2[, c("Longitude", "Latitude")]
 
  CurSpp3 <- SpatialPoints(coordinates, proj4string = crs(df[, "defaultCRS"]))
  
  #Reproject to projection of environmental data
  CurSpp3 <- spTransform(CurSpp3, crs(desiredCRS))
  
  #Extract data
  data <- extract(train2, CurSpp3)
  
  #make it a df and not a vector
  data <- as.data.frame(data)
  
  #add the extracted values back to the CurSpp df
  CurSpp4 <- cbind(Species, CurSpp3@coords, data)
  
  if (length(grep("_", CurSpp)) > 0) {
    s <- as.character(CurSpp4$Species[1])
  } else {
    s <- c(unlist(strsplit(as.character(CurSpp4$Species[1]), " ")))[1]
  }
  
  print(paste0("      Evaluating: ", s))
  
  #write the changes
  s <- gsub(" ", "_" , s)
  #remove any rows with NA's
  CurSpp <- na.omit(CurSpp4)
  
  #write out the .csv with the extracted data for the envSample call
  write.csv(CurSpp, file = paste("samples/", s, ".csv", sep = ""), row.names = FALSE)
  
  if (backgroundPointsStep == "N") {
    
    BPDirList <- list.files(paste0(test, "/backgrounds"), full.names = TRUE)
    CurSpecBPDir <- BPDirList[grep(s, BPDirList)]
    for (f in 1:length(CurSpecBPDir)) {
      BP <- read.csv(file = paste0(CurSpecBPDir[f]))
      if (tolower(names(BP)[1]) == "x" && length(grep("^y$", tolower(names(BP)))) == 0) {
        BP <- BP[, 2:ncol(BP)]
      }
      names(BP)[c(grep("lon", tolower(names(BP))), grep("^x$", tolower(names(BP))))] <- "Longitude"
      names(BP)[c(grep("lat", tolower(names(BP))), grep("^y$", tolower(names(BP))))] <- "Latitude"
      names(BP)[c(grep("^sp", tolower(names(BP))), grep("name", tolower(names(BP))))] <- "Species"
      BP <- BP[, c("Species", "Longitude", "Latitude")]
      BP2 <- BP[, 2:3]
      
      try({
        BPoints <- SpatialPoints(BP2, proj4string = crs(defaultCRS))
        BPoints2 <- spTransform(BPoints, crs(desiredCRS))
      }, silent = TRUE)
      
      if (!exists("BPoints2")) {
        BPoints2 <- SpatialPoints(BPoints, proj4string = crs(desiredCRS))
      }
      
      BPExtract <- extract(train2, BPoints2)
      BP_Final <- cbind(BP[, 1:3], BPExtract)
      BP_Final <- BP_Final[!is.na(BP_Final[, 4]),]
      write.csv(BP_Final, file = paste0(CurSpecBPDir[f]), row.names = FALSE)
      rm(BPoints2)
    }
  }
}

#Run---------------------------------------------------------

setwd(test)
ListSpp <- c()
speciesWorked <- spp_batch
for (i in 1:length(speciesWorked)) {
  ListSpp <- c(ListSpp, list.files(path = "species", pattern = speciesWorked[i], full.names = TRUE))
}
ListSpp <- unique(ListSpp)
print(ListSpp)

nfiles<- length(ListSpp)
dir.create("samples")

BPDirList <- list.files(paste0(test, "/backgrounds"), full.names = TRUE)
if (backgroundPointsStep == "N") {
  for (i in 1:length(ListSpp)) {
    Species1 <- substr(ListSpp[i], 9, nchar(ListSpp[i]) - 4)
    CurSpecBPDir <- BPDirList[grep(Species1, BPDirList)]
    if (length(CurSpecBPDir) == 0) {
      stop(paste0("No background file(s) found for ", Species1))
    }
  }
}

clus <- makeCluster(ncores, outfile = outfile)

clusterExport(clus, varlist = c("ncores", "run", "train", "train2", 
                                "test", "df", "rastertype", "ListSpp", "project", 
                                "defaultCRS", "desiredCRS", "dispersalStep", "backgroundPointsStep"))
clusterEvalQ(clus, library(raster))

print(paste0("   Creating files in:"))
print(paste0("      ", getwd(), "/samples"))
out<-parLapply(clus, ListSpp, function(x) run(x))
stopCluster(clus)
rm(train, train2)
gc()
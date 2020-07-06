####subsampleOccur_bg.R####
##Conducts environmental subsampling on background points (Varela et al. 2014)

#Initializations----------------------------
#Loads the necessary packages
library(dplyr)
library(parallel)
library(raster)
library(rgdal)

#Loads the necessary variables from "df"
test <- df[, "test"]
backgrounds <- "backgrounds"
nsubsamp <- df[, "nsubsamp"]
nclimatebins <- df[, "nclimatebins"]
nPCAxes <- df[, "nPCAxes"]
ncores <- as.numeric(df[, "ncores"])
randomseed <- as.numeric(df[, "ncores"])
outfile <- "statsout.txt"
Categorical <- df[, "Categorical"]

#ptz: df with Latitude, Longitude, extracted environment variables
#no_bins: the number of bins to make in each dimension
#output is a df of latitude, longitude points, one for each occupied bin
set.seed(randomseed)

#Functions-----------------------------------
#Conducts the environmental subsampling	
#ptz: df with Latitude, Longitude, extracted environment variables	
#no_bins: the number of bins to make in each dimension	
#output is a df of latitude, longitude points, one for each occupied bin
envSample <- function (ptz, no_bins = nclimatebins) {
  #Makes a landing spot for the bin membership vectors
  out_ptz <- ptz[, c("Longitude", "Latitude")]
  
  #Runs PCA on ptz if categorical variables are absent
  if (is.na(Categorical)) {
    PCAptz <- prcomp(ptz[, 3:ncol(ptz)], scale = TRUE)
    PCAImp <- summary(PCAptz)$importance
    #Determines the number of PC axes to use for subsampling
    if (!is.na(nPCAxes)) {
      NumberAxes <- nPCAxes
    } else {
      NumberAxes <- max(2, min(which(PCAImp[3,] > 0.95)))
    }
    
    if (NumberAxes > (ncol(ptz) - 3)) {
      NumberAxes <- ncol(ptz) - 3
      message(paste0("'NPCAxes' is more than the number of environmental variables (", ncol(ptz) - 3, "): Setting 'NPCAxes' to ", NumberAxes, "."))
    }
    
    #Adds PCA values to the unsubsampled data frame
    ptz <- cbind(ptz[, 1:2], PCAptz$x[, 1:NumberAxes])
  }
  #Cycles through all of the environmental variable combinations (PC axes) (columns 3 to end)
  for(i in 3:length(names(ptz))) {
    #Makes a data frame that is this variable with no NA values
    k <- ptz[!is.na(ptz[, i]), i]
    #Calculates the observed range of this variable
    rg <- range(k)
    #Figures out the resolution from the number of bins
    res <- (rg[2] - rg[1]) / no_bins
    #Rescales the axis by the range and bin size, so the value is just a 
    #number from 1 to no_bins for its bin membership
    d <- (ptz[, i] - rg[1]) / res
    f <- ceiling(d)
    f[f == 0] <- 1 #Moves the zeros into the 1 bin
    #Corrects the name of the vector, so it will carry over to the output
    names(f) <- names(ptz)[i]
    #Adds the bin membership vector to the output df for this section
    out_ptz <- cbind(out_ptz, f)
    #Gets the names correct
    names(out_ptz)[length(names(out_ptz))] <- names(ptz)[i]
  }
  
  #Subsamples the bin membership df to come up with the filled bins
  sub_ptz <- distinct(out_ptz[, -1:-2])
  #Counts the number of filled bins
  no_grps <- nrow(sub_ptz)
  #Adds a column with the group membership number; this number is arbitrary
  sub_ptz$grp <- c(1:no_grps)
  
  #Joins the out_ptz with the subsample to capture the group membership info
  #note: join() will automatically match the variable names from these two dfs
  out_ptz <- left_join(out_ptz, sub_ptz)
  
  #Selects a random point for each group -- this is an additional improvement on the 
  #Varela et al. function, because it does not pick the same points each time.
  
  #Makes a landing spot for the data
  final_out <- data.frame(x = numeric(), y = numeric())
  
  #Cycles through each group
  for(i in 1:no_grps) {
    
    #Subsets to the members of the ith group, keep only the Latitude and Longitude
    grp_mbrs <- out_ptz[out_ptz$grp == i, c("Longitude", "Latitude")]
    
    #Picks one of these group members to output
    grp_out <- grp_mbrs[sample(1:nrow(grp_mbrs), 1),]
    #Binds this sampled row to the output df
    final_out <- rbind(final_out, grp_out)
  }
  
  #Returns the subsampled points as a df of Latitude and Longitude values
  colnames(final_out)
  return(final_out)
}

run <- function(cur) {
  print("   Writing files to:")
  print(paste0("      ", test, "/", cur))
  
  #Loads the points from the .csv file named 'cur' in the wd and clean (remove NaN and Inf)
  cur2 <- na.omit((read.csv(cur)))
  cur2 <- cur2[!is.infinite(rowSums(cur2[, names(cur2) != "Species"])), ]
  cur2 <- unique(cur2)
  #Pulls the species' name from this file
  if (length(grep("_", cur)) > 0) {
    s <- as.character(cur2$Species[1])
    #Changes 'Species name' to 'Species_name'
    s <- gsub(" ", "_" , s)
  } else {
    s <- c(unlist(strsplit(as.character(cur2$Species[1]), " ")))[1]
  }
  
  #Makes a directory called 'backgrounds/%speciesname%'
  dir.create(paste0("backgrounds/", s))
  
  #Works through the number of subsamples from the config file
  for (r in 1:nsubsamp) {
    #Removes species name column for cleaner subsampling
    coords <- envSample(cur2[, names(cur2) != "Species"])
    #Joins defaults to left join, keeping only rows from the first df
    out <- left_join(coords, cur2) 
    #so 'out' has the structure Longitude, Latitude, Species name, env vars
    #with the rows chosen as the env subsample
    #Reorders columns to be Species, Latitude, Longitude      
    out <- out[, c(3, 1, 2, 4:ncol(out))]
    colnames(out)[2] <- "x"
    colnames(out)[3] <- "y"
    #Writes the subsampled data out
    write.csv(out, file = paste0("backgrounds/", s, "/", s, "_background_", r, ".csv", sep = ""), row.names = FALSE)
  }
} #end of run function.

#Run-------------------------------
#Generates the species list for parallelization
setwd(test)
spp.list <- c()
speciesWorked <- substr(spp_batch, 1, nchar(spp_batch) - 4)
for (i in 1:length(speciesWorked)) {
  spp.list <- c(spp.list, list.files(path = "backgrounds", full.names = TRUE, pattern = paste0(speciesWorked[i], "_background")))
}
spp.list <- unique(spp.list)
print("   Will evaluate species:")
print(spp.list)

#Parallelization
clus <- makeCluster(ncores, outfile = outfile, setup_timeout = 0.5)
clusterExport(clus, varlist = c("randomseed", "test", "samples", "nsubsamp","nclimatebins","spp.list", "nPCAxes", "run", "envSample", "Categorical"))

clusterEvalQ(clus, library(raster))
clusterEvalQ(clus, library(rgdal))
clusterEvalQ(clus, library(dplyr))

out<-parLapply(clus, spp.list, function(x) run(x))
stopCluster(clus)
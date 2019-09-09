#Varela et al. 2014 subsampling technique
#Initializations----------------------------
library(dplyr)
library(parallel)
library(raster)
library(rgdal)

test <- df[, "test"]
samples <- "samples"
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
envSample <- function (ptz, no_bins = nclimatebins) {
  #make a landing spot for the bin membership vectors
  out_ptz <- ptz[, c("Longitude", "Latitude")]
  
  #run PCA on ptz if categorical variables are absent
  if (is.na(Categorical)) {
    PCAptz <- prcomp(ptz[, 3:ncol(ptz)], scale = TRUE)
    PCAImp <- summary(PCAptz)$importance
    #Determine the number of PC axes to use for subsampling
    if (!is.na(nPCAxes)) {
      NumberAxes <- nPCAxes
    } else {
      NumberAxes <- max(2, min(which(PCAImp[3,] > 0.95)))
    }
    
    if (NumberAxes > (ncol(ptz) - 3)) {
      NumberAxes <- ncol(ptz) - 3
      message(paste0("'NPCAxes' is more than the number of environmental variables (", ncol(ptz) - 3, "): Setting 'NPCAxes' to ", NumberAxes, "."))
    }
    
    #Add PCA values to the unsubsampled data frame
    ptz <- cbind(ptz[, 1:2], PCAptz$x[, 1:NumberAxes])
  }
  #cycle through all of the environmental variable combinations (PC axes) (columns 3 to end)
  for(i in 3:length(names(ptz))) {
    #make a data frame that is this variable with no NA values
    k <- ptz[!is.na(ptz[, i]), i]
    #calculate the observed range of this variable
    rg <- range(k)
    #figure out the resolution from the number of bins
    res <- (rg[2] - rg[1]) / no_bins
    #rescale the axis by the range and bin size, so the value is just a 
    #number from 1 to no_bins for its bin membership
    d <- (ptz[, i] - rg[1]) / res
    #d is now a vector of values ranging from 0 to no_bins
    f <- ceiling(d)
    #f is a vector of bin membership
    f[f == 0] <- 1 #move the zeros into the 1 bin
    #correct the name of the vector, so it will carry over to the output
    names(f) <- names(ptz)[i]
    #add the bin membership vector to the output df for this section
    out_ptz <- cbind(out_ptz, f)
    #get the names correct
    names(out_ptz)[length(names(out_ptz))] <- names(ptz)[i]
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
    grp_mbrs <- out_ptz[out_ptz$grp == i, c("Longitude", "Latitude")]
    
    #pick one of these group members to output
    grp_out <- grp_mbrs[sample(1:nrow(grp_mbrs), 1),]
    #bind this sampled row to the output df
    final_out <- rbind(final_out, grp_out)
  }
  
  #return the subsampled points as a df of Latitude and Longitude values
  colnames(final_out)
  return(final_out)
}

run <- function(cur) {
  print("   Writing files to:")
  print(paste0("      ", test, "/", samples, "/", cur))
  
  #load the points from the .csv file named 'cur' in the wd
  cur2 <- na.omit((read.csv(cur)))
  cur2 <- unique(cur2)
  
  #pull the species' name from this file
  if (length(grep("_", cur)) > 0) {
    s <- as.character(cur2$Species[1])
    #change 'Species name' to 'Species_name'
    s <- gsub(" ", "_" , s)
  } else {
    s <- c(unlist(strsplit(as.character(cur2$Species[1]), " ")))[1]
  }
  
  #make a directory called 'samples/%speciesname%'
  dir.create(paste0(samples, "/", s))
  
  #work through the number of subsamples from the config file
  for (r in 1:nsubsamp) {
    #'cur2' comes in with the species name as the first column, so 
    #we'll have to remove that for the envSample() to run.
    coords <- envSample(cur2[, names(cur2) != "Species"])
    #so out is a df of the species name, long and lat of the envSample run
    #with the Longitude called 'x' and the Latitude called 'y'
    #join defaults to left join, keeping only rows from the first df
    out <- left_join(coords, cur2) 
    #so 'out' has the structure Longitude, Latitude, Species name, env vars
    #with the rows chosen as the env subsample
    #need to reorder columns to be Species, Latitude, Longitude      
    out <- out[, c(3, 1, 2, 4:ncol(out))]
    colnames(out)[2] <- "x"
    colnames(out)[3] <- "y"
    #Write the subsampled data out
    write.csv(out,file=paste0(paste0(samples, "/", s, "/"), "OccurrenceSamplePoints_", r, ".csv", sep = ""), row.names = FALSE)
  }
} #end of run function.

#Run-------------------------------
setwd(test)
spp.list <- c()
speciesWorked <- spp_batch
for (i in 1:length(speciesWorked)) {
  spp.list <- c(spp.list, list.files(path = samples, full.names = TRUE, pattern = speciesWorked[i]))
}
spp.list <- unique(spp.list)
print("   Will evaluate species:")
print(spp.list)

clus <- makeCluster(ncores, outfile = outfile)
clusterExport(clus, varlist = c("randomseed", "test", "samples", "nsubsamp","nclimatebins","spp.list", "nPCAxes", "run", "envSample", "Categorical"))

clusterEvalQ(clus, library(raster))
clusterEvalQ(clus, library(rgdal))
clusterEvalQ(clus, library(dplyr))

#here is where it calls the run() function in a parLapply()
out<-parLapply(clus, spp.list, function(x) run(x))
stopCluster(clus)
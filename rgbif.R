####rgbif.R####
##Creates occurrence point files for species by interfacing with GBIF

#Initializations-----------------------------------
#Loads the necessary packages
library(rgbif)

#Loads the necessary variables from "df"
spplist <- df[, "spplist"]
occurrences <- df[, "occurrences"]
sppcountsloc <- df[, "sppcountsloc"]
lat <- df[, "decimalLatitude"] 
long <- df[, "decimalLongitude"]
minlat <- as.numeric(df[, "minlat"])
maxlat <- as.numeric(df[, "maxlat"])
minlong <- as.numeric(df[, "minlong"])
maxlong <- as.numeric(df[, "maxlong"])
test <- df[, "test"]
result_dir <- df[, "result_dir"]
dispersalStep <- df[, "dispersalStep"]
if (dispersalStep == "Y") {
  dispersalRate_dir <- df[, "dispersalRate_dir"]
}
#Allows for connection to last longer than 60 seconds
options(timeout = 1000)
#Copies the maxent javascript to the test folder
setwd(test)
file.copy(paste0(DataDirectory, "/maxent.jar"), test)

#Reads species list csv file
SppList <- read.csv(paste0(spplist), strip.white = TRUE, stringsAsFactors = FALSE)
OurSpp <- SppList
OurSpp <- data.frame(lapply(OurSpp, as.character),stringsAsFactors = FALSE)

#Standardizes headings
names(OurSpp)[1] <- "Taxon"
names(OurSpp)[2] <- "Scientific.Name"

#Creates new data frame with taxonomy, number of occurrences, etc.
nspp <- nrow(OurSpp)
OurSpp <- as.data.frame(cbind(SpeciesSearched = rep(NA, times = nspp), 
                              OurSpp, 
                              OrigOccurrences = rep(0, times = nspp), 
                              Occurrences = rep(0, times = nspp), 
                              StudyArea_Occur = rep(0, nspp), 
                              Keys = rep(NA, times = nspp), 
                              Family = rep(NA, times = nspp), 
                              Genus = rep(NA, times = nspp)))

setwd(occurrences)

#Taxonomy Matching and Key Generation-----------------------------------------
#Matches species to GBIF taxonomy and determines if names need to be changed and/or merged
for (f in 1:nspp) {
  
  #Gets the GBIF backbone for each scientific name
  s <- as.character(OurSpp$Scientific.Name[f])
  SpeciesName <- name_backbone(s, strict = TRUE)
  
  #If the scientific name is not in the GBIF backbone, prints an error
  if (SpeciesName$matchType == "NONE") {
    message("Species ", s, " not found: Check spelling or try searching a synonym")
    next
  }
  
  SpecSplit <- unlist(strsplit(s, " "))
  
  #If the taxon is a species (not a subspecies), renames it to the GBIF taxonomy.
  #However, if the taxon *is* a subspecies, leave alone.
  if (length(SpecSplit) == 2) {
    if (SpeciesName$species != s) {
      
      #If the GBIF backbone name doesn't equal the species name provided, either renames or merges with another provided species
      if (length(grep(SpeciesName$species, OurSpp$Scientific.Name) > 0)) {
        message(paste0(s, " is the same species as ", SpeciesName$species, " (also listed in spplist), and will therefore be discarded from this analysis"))
      } else {
        message(paste0(s, " will be renamed to: ", SpeciesName$species))
      }
      
      s <- SpeciesName$species
      
      #Renames species dispersal rate data (if applicable)
      if (dispersalStep == "Y") {
        #Loads dispersal rate data and find row with matching species name
        DispDataName <- list.files(path = dispersalRate_dir, pattern = ".csv")
        DispData <- read.csv(paste0(dispersalRate_dir, "/", DispDataName[1]), stringsAsFactors = FALSE)
        DispSpec <- grep(" ", DispData[1, ])
        Specloc <- grep(paste0(s, "$"), DispData[, DispSpec])
        DispData[Specloc, DispSpec] <- s
      }
      
      #Renames species lists and adds key to "OurSpp" dataframe
      OurSpp$Scientific.Name[f] <- s
      SppList[f, 2] <- s
    }
    #Adds species key to "OurSpp" dataframe
    OurSpp$Keys[f] <- SpeciesName$speciesKey
  } else {
    #If there is an "acceptedUsageKey" in the dataframe, renames to that subspecies name. 
    #Otherwise, leaves it alone, with a warning if the given species name and the GBIF taxonomy don't match
    if ("acceptedUsageKey" %in% colnames(SpeciesName)) {
      
      #Find the preferred subspecific name for the given species and renames it
      ss_search <- occ_search(SpeciesName$acceptedUsageKey, limit = 1, fields = c('species', 'infraspecificEpithet'))
      new_ss <- paste0(ss_search$data$species, " ", ss_search$data$infraspecificEpithet)
      message(paste0(s, " will be renamed to: ", new_ss))
      
      s <- SpeciesName$species
      
      #Renames species dispersal rate data (if applicable)
      if (dispersalStep == "Y") {
        #Loads dispersal rate data and find row with matching species name
        DispDataName <- list.files(path = dispersalRate_dir, pattern = ".csv")
        DispData <- read.csv(paste0(dispersalRate_dir, "/", DispDataName[1]), stringsAsFactors = FALSE)
        DispSpec <- grep(" ", DispData[1, ])
        Specloc <- grep(paste0(s, "$"), DispData[, DispSpec])
        DispData[Specloc, DispSpec] <- s
      }
      
      #Renames species lists
      OurSpp$Scientific.Name[f] <- s
      SppList[f, 2] <- s
      #Adds species key to "OurSpp" dataframe
      OurSpp$Keys[f] <- SpeciesName$acceptedUsageKey
    } else {
      if (SpeciesName$species != paste(SpecSplit[1:2], collapse = " ")) {
        message(paste0("Warning: '", s, "' has a different species name than the accepted GBIF taxonomy: ", SpeciesName$species))
        message("         If fewer records than expected are returned, check the GBIF taxonomy for the preferred suspecies name")
      }   
      #Adds species key to "OurSpp" dataframe
      OurSpp$Keys[f] <- SpeciesName$usageKey
    }
  }
}

#Creates species list for GBIF scraping
#This also removes any duplicates that occur because of the taxonomy changes
SppList <- unique(SppList)
OurSpp <- unique(OurSpp)
OurSpp <- OurSpp[which(SppList[, 2] != "NA"), ]
SppList <- SppList[which(SppList[, 2] != "NA"), ]
nspp <- nrow(OurSpp)

#Extracting Occurrence Data----------------------------------------
#Suppresses URL warnings
options(warn = -1)

#Extracts occurrences from GBIF

FailedSpecies <- c()
#Iterates through all taxa
speciterate <- function(start, finish) {
  for(i in start:finish) {
    p <<- i
    #Gets taxon name
    s <- as.character(OurSpp$Scientific.Name[i]) 
    SpecSplit <- unlist(strsplit(s, " "))
    
    #Determines if the taxon is a species or a subspecies and splits workflow accordingly 
    if (length(SpecSplit) == 2) {
      #Prints progress report
      print(paste("Species ", i, " of ", nspp, ": ", s, sep = ""))
      print(paste0("   Beginning search: ", Sys.time()))
      if (is.na(OurSpp$Keys[i])) {
        message(paste0("Species ", s, " not found: check spelling or try a synonym"))
        FailedSpecies <- c(FailedSpecies, s)
        next()
      }
      OurSpp$OrigOccurrences[i] <- 0
      OurSpp$Occurrences[i] <- 0
      
      #Searches for occurrences points given the species keys
      Occ <- occ_search(taxonKey = OurSpp$Keys[i], 
                        decimalLatitude = lat, 
                        decimalLongitude = long, 
                        hasCoordinate = TRUE, 
                        limit = 200000, 
                        fields = c('species', 
                                   'infraspecificEpithet', 'decimalLatitude', 
                                   'decimalLongitude', 'basisOfRecord', 
                                   'issues','locality', 'elevation', 
                                   'elevationAccuracy', 'continent', 
                                   'stateProvince', 'county',
                                   'year', 'month', 'day', 'evenDate', 
                                   'references', 'license', 'geodeticDatum', 
                                   'gbifID', 'type', 'preparations', 
                                   'catalogNumber', 'occurrenceStatus'))$data
      
      if(is.null(Occ)) {
        message("No occurrences found within study area! Check species name or study area extent")
        FailedSpecies <- c(FailedSpecies, s)
        next()
      }
      
      OurSpp$OrigOccurrences[i] <- nrow(Occ)
      
    } else {
      print(paste("Subspecies ", i, " of ", nspp, ": ", s, sep = ""))
      print(paste0("   Beginning search: ", Sys.time()))
      OurSpp$OrigOccurrences[i] <- 0
      OurSpp$Occurrences[i] <- 0
      
      #Searches for occurrences points given the species keys
      Occ <- occ_search(taxonKey = OurSpp$Keys[i], 
                        decimalLatitude = lat, 
                        decimalLongitude = long, 
                        hasCoordinate = TRUE, 
                        limit = 200000, 
                        fields = c('species', 
                                   'infraspecificEpithet', 'decimalLatitude', 
                                   'decimalLongitude', 'basisOfRecord', 
                                   'issues','locality', 'elevation', 
                                   'elevationAccuracy', 'continent', 
                                   'stateProvince', 'county',
                                   'year', 'month', 'day', 'evenDate', 
                                   'references', 'license', 'geodeticDatum', 
                                   'gbifID', 'type', 'preparations', 
                                   'catalogNumber', 'occurrenceStatus'))$data
      
      
      if(is.null(Occ)) {
        message("No occurrences found within study area! Check species name or study area extent")
        FailedSpecies <- c(FailedSpecies, s)
        next()
      }
      
      #Adds subspecies name to the species column and deletes subspecies column
      Occ$species <- paste0(Occ$species, " ", Occ$infraspecificEpithet)
      OurSpp$OrigOccurrences[i] <- nrow(Occ)
    }
    
    #Print progress to console
    print(paste0("   Finishing search: ", Sys.time()))
    print(paste0("   Number original occurrences: ", nrow(Occ)))
    
    if(nrow(Occ) == 0) {
      print(paste0("   Species failed, no search data found: ", s))
    } else if (!is.atomic(Occ)) {
      #Removes fossils
      Occ <- Occ[!Occ$basisOfRecord == "FOSSIL_SPECIMEN", ]
      #Removes geographical issues with the data
      Occ <- Occ[!grepl("cdiv",  Occ$issues), ]
      Occ <- Occ[!grepl("cdout",  Occ$issues), ]
      Occ <- Occ[!grepl("cdrepf",  Occ$issues), ]
      Occ <- Occ[!grepl("cdreps",  Occ$issues), ]
      Occ <- Occ[!grepl("gdativ",  Occ$issues), ]
      Occ <- Occ[!grepl("preneglat",  Occ$issues), ]
      Occ <- Occ[!grepl("preneglon",  Occ$issues), ]
      Occ <- Occ[!grepl("preswcd",  Occ$issues), ]
      Occ <- Occ[!grepl("txmatnon",  Occ$issues), ]
      Occ <- Occ[!grepl("zerocd",  Occ$issues), ]
      #Removes duplicates, log those removed
      Occ <- Occ[!duplicated(data.frame(Occ$decimalLatitude, Occ$decimalLongitude)), ]
      OurSpp$Occurrences[i] <- nrow(Occ)
      
      #clips the occurrences to the study area, logs species counts
      SA_Occ <- Occ[as.numeric(Occ$decimalLatitude) <= maxlat, ]
      SA_Occ <-SA_Occ[as.numeric(SA_Occ$decimalLatitude) >= minlat, ]
      SA_Occ <-SA_Occ[as.numeric(SA_Occ$decimalLongitude) <= maxlong, ]
      SA_Occ<- SA_Occ[as.numeric(SA_Occ$decimalLongitude) >= minlong, ]
      OurSpp$StudyArea_Occur[i] <- nrow(SA_Occ) 
      print(paste0("   Number in the study area: ", nrow(SA_Occ)))
      OurSpp$SpeciesSearched[i] <- s
      
      #Records details about species as found by GBIF
      sppkeys <- OurSpp$Keys[i]
      specurl <- paste("http://api.gbif.org/v1/species/", trimws(sppkeys[1]), sep="")
      URLRead <- colnames(read.csv(specurl))
      gbifapidata <- gsub("\\.", "_", URLRead)
      OurSpp$Species[i] <- paste(as.character(substr(gbifapidata[grep("^species_", gbifapidata)], 9, nchar(gbifapidata[grep("^species_", gbifapidata)]))))
      
      if (length(grep("^class_", gbifapidata)) > 0) {
        OurSpp$Class[i] <- paste(as.character(substr(gbifapidata[grep("^class_", gbifapidata)], 7, nchar(gbifapidata[grep("^class_", gbifapidata)]) - 1)))
      }
      if (length(grep("^family_", gbifapidata)) > 0) {
        OurSpp$Family[i] <- paste(as.character(substr(gbifapidata[grep("^family_", gbifapidata)], 8, nchar(gbifapidata[grep("^family_", gbifapidata)]))))
      }
      if (length(grep("^genus_", gbifapidata)) > 0) {
        OurSpp$Genus[i] <- paste(as.character(substr(gbifapidata[grep("^genus_", gbifapidata)], 7, nchar(gbifapidata[grep("^genus_", gbifapidata)]))))
      }
      
      #Writes occurrences
      setwd(occurrences)
      #Uses species name from our SppList file as occurrence file name for the csv
      write.csv(Occ, file = paste(gsub(" ", "_", s), ".csv", sep = "")) 
      print(paste0("   Finishing species: ", Sys.time()))
    } else {
      print(paste0("   Species failed, no search data found: ", s))
    }
    
    OurSpp <<- OurSpp
    #Creates csv file with list of species
    #Adds species with 0 occurrences to failed species
    if (OurSpp$Occurrences[i] == 0) {
      FailedSpecies <- c(FailedSpecies, s)
    }
    FailedSpecies <<- FailedSpecies
    rm(Occ)
    gc()
  }
}

#It is looped to avoid timeout or temporary internet connectivity issues
speciterate(1, nspp)
if (exists("p")) {
  while (p < nspp) {
    speciterate(p, nspp)
  }
}

#File copying and management--------------------
convert <- function(CurSpp) {
  tryCatch({
    #Ensures that species names are consistent with subsequent steps
    CurSpp2 <- (read.csv(CurSpp))
    CurSpp2 <- data.frame(lapply(CurSpp2, as.character), stringsAsFactors = FALSE)
    s <- as.character(CurSpp2$species[1])
    
    #If there is more than one taxon name (usually due to taxonomic discrepancies), make everything the same name
    if (length(unique(CurSpp2$species)) > 1) {
      CurSpp2$species <- s
    }
    
    #Renames column headings to be used in subsequent steps
    CurSppFinal <- CurSpp2[which(CurSpp2$species == s), ]
    names(CurSppFinal)[names(CurSppFinal) == "species"] <- "Species"
    names(CurSppFinal)[names(CurSppFinal) == "decimalLongitude"] <- "Longitude"
    names(CurSppFinal)[names(CurSppFinal) == "decimalLatitude"] <- "Latitude"
    CurSppFinal <- CurSppFinal[, c("Species", "Longitude", "Latitude")]
    
    sSpace <- gsub(" ", "_", s)
    write.csv(CurSppFinal, file = paste0("species/", sSpace, ".csv"), row.names = FALSE)
  }, error = function(err) {
    print(paste("MY_ERROR: ", CurSpp, " ", err))
    return(paste0("error: ", err))
  })
}

#Write the species counts data frame to CSV file
setwd(sppcountsloc)
write.csv(OurSpp, file = df[, "counts"])

#Copies occurrence files to the occurrences folder
setwd(occurrences)
ListSpp <- list.files(path = getwd(), pattern = '\\.csv', full.names = TRUE)
nspp <- length(ListSpp)
for(j in 1:nspp) {
  file.copy(paste0(ListSpp), test)
}

#Create species directory
setwd(test)
if (!dir.exists(paste0(test, "/species"))) {
  dir.create("species")
}
ListSpp <- list.files(path = getwd(), pattern = '\\.csv', full.names = TRUE)
ListSpp <- ListSpp[1:length(ListSpp)]

#Convert and write occurrence data to be used in megaSDM
out<- lapply(ListSpp, function(x) convert(x))
write.csv(SppList, file = paste0(df[, "spplist"]), row.names = FALSE)

#Write "FailedSpecies" csv
if (length(FailedSpecies) > 0) {
  message(paste0("Species Generated in Failed Species List: Check ", paste0(result_dir, "/", "FailedSpecies.csv"), " for failed species"))
}
write.csv(FailedSpecies, paste0(result_dir, "/", "FailedSpecies.csv"))

#Re-writes out changes dispersal file (if species names needed to be changed)
if (!exists("DispData") && dispersalStep == "Y") {
  DispDataName <- list.files(path = dispersalRate_dir, pattern = ".csv")
  DispData <- read.csv(paste0(dispersalRate_dir, "/", DispDataName[1]), stringsAsFactors = FALSE)
}
if (dispersalStep == "Y") {
  write.csv(DispData, file = paste0(dispersalRate_dir, "/", DispDataName), row.names = FALSE)
}
gc()
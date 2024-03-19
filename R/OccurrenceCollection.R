#' Download and vet GBIF occurrence data
#'
#' Takes a list of species and collects occurrence data from GBIF (Global Biodiversity Information Facility).
#' Acts as a wrapper for \code{rgbif::occ_search}; however, this function is much more efficient for a large
#' number of species. It also checks the taxonomy of the given species list against the GBIF taxonomy,
#' renaming or merging taxa if necessary. Furthermore, this function vets the occurrence data, removing
#' occurrence points that are of insufficient quality for species distribution modelling.
#' Finally, \code{OccurrenceCollection()}
#' provides the number of occurrences found within given training and study areas.
#'
#' @param spplist a vector of scientific names, using GBIF taxonomy. Names can be species or subspecies.
#' @param output A full directory name where the downloaded species occurrences will be written to.
#' @param trainingarea Extent object, or vector of desired training extent in form \code{c(xmin, xmax, ymin, ymax)}.
#' Given in latlong coordinates. If set to "NA" (the default), all occurrence points will be
#' downloaded (up to 99,999 points), regardless of their location.
#' @param studyarea (optional) Extent object, or vector of the area that the SDM will be projected on in
#' form \code{c(xmin, xmax,ymin, ymax)}. Given in latlong coordinates. If provided, the number of occurrence
#' points found within this study region will be calculated.
#' @export
#' @return Writes .csv files of GBIF occurrences to a directory provided by \code{output}. If any species failed
#' (for example, the scientific name was not found in the search or no occurrences exist within a provided \code{trainingarea}),
#' another .csv file is written out in the same folder with the names of the species that failed. In addition,
#' a dataframe is returned by the function that contains the high taxonomy of each species and the nubmer of
#' occurrences found within the \code{trainingarea} and, optionally, the \code{studyarea}.
#'
#' List of occurrences removed by \code{OccurrenceCollection}:
#' \enumerate{
#' \item Fossil specimens (mitigates the effect of long-term climiatic changes on the SDM)
#' \item "cdiv": Coordinate Invalid
#' \item "cdout": Coordinate Out of Range
#' \item "cdrepf": Coordinate Reprojection Failed
#' \item "cdreps": Coordinate Reprojection Suspicious
#' \item "gdativ": Geodetic Datum Invalid
#' \item "preneglat": Presumed Negated Latitude
#' \item "preneglon": Presumed Negated Longitude
#' \item "preswcd": Presumed Swapped Coordinates
#' \item "txmatnon": No Taxon Match
#' \item "zeocd": Exact 0/0 Coordinate.}
#' Further vetting may be done by hand.

OccurrenceCollection <- function(spplist,
                                 output,
                                 trainingarea = NA,
                                 studyarea = NA) {
  options(timeout = 1000)

  if (!dir.exists(output)) {
    dir.create(output)
  }

  #Creates new data frame with taxonomy, number of occurrences, etc.
  nspp <- length(spplist)
  OurSpp <- data.frame(SpeciesSearched = rep(NA, times = nspp),
                                Scientific.Name = as.character(sub("_", " ", spplist)),
                                OrigOccurrences = rep(0, times = nspp),
                                Occurrences = rep(0, times = nspp),
                                StudyArea_Occur = rep(0, nspp),
                                Keys = as.numeric(rep(NA, times = nspp)),
                                Class = rep(NA, times = nspp),
                                Family = rep(NA, times = nspp),
                                Genus = rep(NA, times = nspp),
                                stringsAsFactors = FALSE)

  #Matches species to GBIF taxonomy and determines if names need to be changed and/or merged
  for (f in 1:nspp) {

    #Gets the GBIF backbone for each scientific name
    s <- as.character(OurSpp$Scientific.Name[f])
    SpeciesName <- rgbif::name_backbone(s, strict = TRUE)

    #If the scientific name is not in the GBIF backbone, prints an error
    if (SpeciesName$matchType == "NONE") {
      message("Species ", s, " not found: Check spelling or try searching a synonym")
      next
    }

    SpecSplit <- unlist(strsplit(s, " "))

    #If the taxon is a species (not a subspecies), renames it to the GBIF taxonomy.
    #However, if the taxon *is* a subspecies, leave alone.
    if (length(SpecSplit) == 2) {
      if(!("species" %in% names(SpeciesName))) {
        message("Species ", s, " not found: Check spelling or try searching a synonym")
        next
      }
      if (SpeciesName$species != s) {

        #If the GBIF backbone name doesn't equal the species name provided, either renames or merges with another provided species
        if (length(grep(SpeciesName$species, OurSpp$Scientific.Name) > 0)) {
          message(paste0(s, " is the same species as ", SpeciesName$species, " (also listed in spplist), and will therefore be discarded from this analysis"))
        } else {
          message(paste0(s, " will be renamed to: ", SpeciesName$species))
        }

        s <- SpeciesName$species

        #Renames species lists and adds key to "OurSpp" dataframe
        OurSpp$Scientific.Name[f] <- s
        spplist[f] <- s
      }
      #Adds species key to "OurSpp" dataframe
      OurSpp$Keys[f] <- SpeciesName$speciesKey
    } else {
      #If there is an "acceptedUsageKey" in the dataframe, renames to that subspecies name.
      #Otherwise, leaves it alone, with a warning if the given species name and the GBIF taxonomy don't match
      if ("acceptedUsageKey" %in% colnames(SpeciesName)) {

        #Find the preferred subspecific name for the given species and renames it
        ss_search <- rgbif::occ_search(SpeciesName$acceptedUsageKey, limit = 1, fields = c('species', 'infraspecificEpithet'))
        new_ss <- paste0(ss_search$data$species, " ", ss_search$data$infraspecificEpithet)
        message(paste0(s, " will be renamed to: ", new_ss))

        s <- SpeciesName$species

        #Renames species lists
        OurSpp$Scientific.Name[f] <- new_ss
        spplist[f] <- new_ss
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

  #Removes any duplicates that occur because of the taxonomy changes
  OurSpp <- unique(OurSpp)
  nspp <- nrow(OurSpp)

  if (!is.na(trainingarea)[1]){
    #Uploads training data and gets extent of occurrence sampling
    decimalLongitude <- c(trainingarea[1:2])
    decimalLatitude <- c(trainingarea[3:4])
  } else {
    decimalLongitude <- c(-180, 180)
    decimalLatitude <- c(-90, 90)
  }

  if (!is.na(studyarea)) {
    studyLongitude <- c(studyarea[1:2])
    studyLatitude <- c(studyarea[3:4])
  }

  #Suppresses URL warnings and sets a vector for failed species
  options(warn = -1)
  FailedSpecies <- c()

  #Iterates through all taxa
  speciterate <- function(start, finish) {
    p <- c()
    for(i in start:finish) {
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
        Occ <- rgbif::occ_search(taxonKey = OurSpp$Keys[i],
                          decimalLatitude = paste(decimalLatitude, collapse = ","),
                          decimalLongitude = paste(decimalLongitude, collapse = ","),
                          hasCoordinate = TRUE,
                          limit = 99999,
                          fields = c('species',
                                     'decimalLatitude',
                                     'decimalLongitude', 'basisOfRecord',
                                     'issues','locality', 'elevation',
                                     'elevationAccuracy', 'continent',
                                     'stateProvince', 'county',
                                     'year', 'month', 'day', 'evenDate',
                                     'references', 'license', 'geodeticDatum',
                                     'gbifID', 'type', 'preparations',
                                     'catalogNumber', 'occurrenceStatus', 'datasetKey'))$data
        # Got rid of 'infraspecificEpithet' because it only exists for subspecies
        if(is.null(Occ)) {
          message("No occurrences found within study area! Check species name or study area extent")
          FailedSpecies <- c(FailedSpecies, s)
          p <- i
          next()
        }

        OurSpp$OrigOccurrences[i] <- nrow(Occ)

      } else {
        print(paste("Subspecies ", i, " of ", nspp, ": ", s, sep = ""))
        print(paste0("   Beginning search: ", Sys.time()))
        OurSpp$OrigOccurrences[i] <- 0
        OurSpp$Occurrences[i] <- 0

        #Searches for occurrences points given the species keys
        Occ <- rgbif::occ_search(taxonKey = OurSpp$Keys[i],
                          decimalLatitude = paste(decimalLatitude, collapse = ","),
                          decimalLongitude = paste(decimalLongitude, collapse = ","),
                          hasCoordinate = TRUE,
                          limit = 99999,
                          fields = c('species',
                                     'infraspecificEpithet', 'decimalLatitude',
                                     'decimalLongitude', 'basisOfRecord',
                                     'issues','locality', 'elevation',
                                     'elevationAccuracy', 'continent',
                                     'stateProvince', 'county',
                                     'year', 'month', 'day', 'evenDate',
                                     'references', 'license', 'geodeticDatum',
                                     'gbifID', 'type', 'preparations',
                                     'catalogNumber', 'occurrenceStatus', 'datasetKey'))$data


        if(is.null(Occ)) {
          message("No occurrences found within study area! Check species name or study area extent")
          FailedSpecies <- c(FailedSpecies, s)
          p <- i
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

        if (exists("studyLongitude")) {
          #clips the occurrences to the study area, logs species counts
          SA_Occ <- Occ[as.numeric(Occ$decimalLongitude) >= studyLongitude[1], ]
          SA_Occ <- SA_Occ[as.numeric(SA_Occ$decimalLongitude) <= studyLongitude[2], ]
          SA_Occ <- SA_Occ[as.numeric(SA_Occ$decimalLatitude) >= studyLatitude[1], ]
          SA_Occ <- SA_Occ[as.numeric(SA_Occ$decimalLatitude) <= studyLatitude[2], ]
          OurSpp$StudyArea_Occur[i] <- nrow(SA_Occ)
          print(paste0("   Number in the study area: ", nrow(SA_Occ)))
        } else {
          OurSpp$StudyArea_Occur[i] <- NA
        }

        OurSpp$SpeciesSearched[i] <- s

        #Records details about species as found by GBIF
        sppkeys <- OurSpp$Keys[i]
        specurl <- paste("http://api.gbif.org/v1/species/", trimws(sppkeys[1]), sep="")
        URLRead <- colnames(utils::read.csv(specurl))
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


        #Uses species name from our spplist file as occurrence file name for the csv
        utils::write.csv(Occ, file = file.path(output, paste0(gsub(" ", "_", s), ".csv")), row.names = FALSE)
        print(paste0("   Finishing species: ", Sys.time()))
      } else {
        print(paste0("   Species failed, no search data found: ", s))
      }

      OurSpp <<- OurSpp
      FailedSpecies <<- FailedSpecies
      rm(Occ)
      gc()
      p <- i
    }
    return(p)
  }

  #It is looped to avoid timeout or temporary internet connectivity issues
  TestNum <- speciterate(1, nspp)
  if (is.null(TestNum)) {
    TestNum <- nspp
  }

  if (TestNum < nspp) {
    while (TestNum < nspp) {
      TestNum <- speciterate(TestNum, nspp)
    }
  }

  #Adds species with 0 occurrences to failed species
  FailedSpecies <- unique(c(FailedSpecies, OurSpp$Scientific.Name[which(OurSpp$OrigOccurrences == 0)]))

  #Write "FailedSpecies" csv
  if (length(FailedSpecies) > 0) {
    message(paste0("Species Generated in Failed Species List: Check ", paste0(output, "/", "FailedSpecies.csv"), " for failed species"))
    utils::write.csv(FailedSpecies, paste0(output, "/", "FailedSpecies.csv"), row.names = FALSE)
  }
  OurSpp$OrigOccurrences <- as.numeric(OurSpp$OrigOccurrences)
  OurSpp$Occurrences <- as.numeric(OurSpp$Occurrences)
  OurSpp$StudyArea_Occur <- as.numeric(OurSpp$StudyArea_Occur)

  if (!dir.exists(paste0(output, "/SpeciesCounts"))) {
    dir.create(paste0(output, "/SpeciesCounts"))
  }

  utils::write.csv(OurSpp, paste0(output, "/SpeciesCounts/SpeciesCounts.csv"), row.names = FALSE)
  return(OurSpp)
}

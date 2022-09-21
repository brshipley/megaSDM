#' Use species-specific sets of environmental data for SDMs
#'
#' Using different sets of environmental variables to model species
#' distributions can help to make more informative species distribution
#' models. This function allows for the modelling of each species based
#' upon a unique subset of the environmental variables.
#'
#' @param occlist a list of .csv files, each containing the occurrence
#' points of a given taxon. The files should be named the same as the taxon
#' of interest (e.g.,: ".../Canis_lupus.csv").
#' @param bglist a list of .csv files corresponding to the background points
#' provided for each taxon. This list of files should be order in the same
#' way as the files given by \code{occlist}.
#' @param env_vars a vector containing the names of the environmental variables
#' (must exactly match column names) to keep for each species as a character
#' strings separated by commas. For example, \code{c("Bio1,Bio3", "Bio1,Bio12")}.
#' Like above, this vector should be in the same order as the lists of
#' species/background points.
#' @param occ_output the directory where output occurrence files will be placed.
#' @param bg_output the directory where output background files will be placed.
#' @export
#' @return Writes .csv files of species occurrences to the sub-directory given by the \code{occ_output} argument
#' and .csv files of background points to the sub-directory given by the \code{bg_output} argument.

VariableEnv <- function(occlist, bglist, env_vars,
                        occ_output, bg_output) {

  #Create sub-directories for background points and occurrence points
  if (!dir.exists(occ_output)) {
    dir.create(occ_output)
  }

  if (!dir.exists(bg_output)) {
    dir.create(bg_output)
  }

  for (i in 1:length(occlist)) {
    CurSpp <- occlist[i]
    s <- grep(CurSpp, occlist)

    #Determines the species name from the name of the occurrence file
    SpeciesSplit <- unlist(strsplit(occlist[s], "/", fixed = TRUE))
    SpeciesName <- substr(SpeciesSplit[length(SpeciesSplit)], 1,
                          nchar(SpeciesSplit[length(SpeciesSplit)]) - 4)

    OccurrenceFile <- occlist[s]
    BackgroundFile <- bglist[s]

    #Determines which environmental variables were requested for that species
    ReqVarList <- as.vector(env_vars[s])
    ReqVarList <- unlist(strsplit(ReqVarList, ",\\s*", perl = TRUE))
    ReqVar <- paste(paste0("^", ReqVarList, "$"), collapse = "|")

    #Subset the background files to only include the required environmental variables
    BGFile <- utils::read.csv(BackgroundFile)
    BGFileCoords <- BGFile[, c("Species", "x", "y")]
    ReqVarCol <- grep(tolower(ReqVar), tolower(names(BGFile)))
    if (length(ReqVarCol) != length(ReqVarList)) {
      stop(paste0("One or more climate variables given in 'env_vars' do not correspond to provided environmental rasters", "\n",
                  "          Check the environmental variable list in 'env_vars'"))
    }

    BGFile_Final <- data.frame(BGFileCoords, BGFile[, ReqVarCol])
    utils::write.csv(BGFile_Final, paste0(bg_output, "/", SpeciesName, "_background.csv"), row.names = FALSE)

    #Subset the occurrence files to only include the required environmental variables
    OCFile <- utils::read.csv(OccurrenceFile)
    OCFileCoords <- OCFile[, c("Species", "x", "y")]
    ReqVarCol <- grep(tolower(ReqVar), tolower(names(OCFile)))
    if (length(ReqVarCol) != length(ReqVarList)) {
      stop(paste0("One or more climate variables given in 'env_vars' do not correspond to provided environmental rasters", "\n",
                  "          Check the environmental variable list in 'env_vars'"))
    }

    OCFile_Final <- data.frame(OCFileCoords, OCFile[, ReqVarCol])
    utils::write.csv(OCFile_Final, paste0(occ_output, "/", SpeciesName, ".csv"), row.names = FALSE)

  }
}

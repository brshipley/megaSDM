#' Generate null distribution models for effective AUC comparison
#'
#' AUC values rely on both omission rate (false negatives) and commision rate (false positives);
#' however, MaxEnt is a presence-only method, making raw AUC values uninformative for
#' comparing across models and species (Jimenez-Valverde 2012). One way to use AUC values to
#' examine presence-only model predictions is to generate model replicates using randomly-generated
#' occurrence data, evaluating their performance using a subset of the real occurrence data.
#' This function generates null models and calculates the Test AUC values when applied to
#' the subset of real occurrence data, for comparison with the model training on the actual
#' data. This method was developed by Bohl et al. 2019.
#'
#' @param envdata a SpatRaster or list of raster files corresponding to the
#' area the model will be trained on.
#' @param replicates how many times should the null model be run to get a distribution
#' of AUC values? Default is 50 replicates.
#' @param bufflist (optional) if background points were spatially-constrained, provide
#' the paths to the buffer files (.shp) used.
#' @param modelpar a list of the arguments passed from the \code{MaxEntModel} function.
#' These arguments should be exactly the same as the actual model generation for effective
#' comparison. For the method described by Bohl et al. (2019), a list of test-samples must
#' be given (evaluate models on the real data). Otherwise, the null distribution will be
#' evaluated on a subset of the random samples (see Raes & ter Steege (2007))
#' @export
#' @return A .csv file (NullModel_AUC.csv) with the Test AUC values for each replicate
#' of the null model.

nullAUC <- function(envdata, replicates = 50, bufflist = NA, modelpar) {

  NullMaxent <- function(occlist, bglist, model_output,
                          ncores = 1, nrep = 1, categorical = NA,
                          alloutputs = TRUE, reptype = "Subsample",
                          test_percent = 20, hinge = TRUE,
                          testsamples = FALSE, regularization = 1) {
    if (length(occlist) < ncores) {
      ncores <- length(occlist)
    }
    ListSpp <- matrix(data = occlist, ncol = ncores)

    hinge <- tolower(hinge)
    alloutputs <- tolower(alloutputs)

    checkError <- function(result, spp.name, failed_runs, run) {
      if (grepl("error", result)) {
        cur <- paste0(spp.name, ", failed run: ", run)
        failed_runs <- c(failed_runs, cur)
        return(failed_runs)
      } else {
        return(failed_runs)
      }
    }

    if(class(envdata) != "SpatRaster") {
      envdata <- terra::rast(envdata)
    }
    
    run <- function(CurSpp) {
      s <- grep(CurSpp, occlist)

      OccurrenceFile <- occlist[s]
      BackgroundFile <- bglist[s]
      TestFile <- testsamples[s]
      
      #Determines the species name from the name of the occurrence file
      SpeciesSplit <- unlist(strsplit(OccurrenceFile, "/", fixed = TRUE))
      spp.name <- substr(SpeciesSplit[length(SpeciesSplit)], 1,
                         nchar(SpeciesSplit[length(SpeciesSplit)]) - 4)

      #Creates sub-directory for the given species
      if (!dir.exists(file.path(model_output, spp.name))) {
        dir.create(file.path(model_output, spp.name))
      }

      failed_runs <- c()

      dir.create(file.path(model_output, spp.name, "NullModels"))
      NullAUCVals <- c()
      if (!is.na(bglist[s])) {
        BGPoints <- read.csv(bglist[s])
        BufferFile <- terra::vect(BGPoints, geom=c("x", "y"))
      }
      NOcc <- nrow(utils::read.csv(OccurrenceFile))
      for (i in 1:replicates) {
        if(!is.na(bglist[s])) {
          RandomSP <- sample(BufferFile, NOcc)
          SPCoords <- data.frame(terra::geom(RandomSP))
          RandomOcc <- terra::extract(envdata, RandomSP)
          RandomOcc <- data.frame(SPCoords[, c("x", "y")], RandomOcc)
          RandomOcc <- data.frame(RandomOcc[stats::complete.cases(RandomOcc),])
        } else {
          RandomOcc <- terra::spatSample(envdata, size = NOcc, method = "random", na.rm = TRUE, xy = TRUE)
        }
        RandomOcc <- data.frame(Species = rep(spp.name, nrow(RandomOcc)), RandomOcc)
        BGFile <- utils::read.csv(BackgroundFile)
        BCol <- c()
        for (b in 1:ncol(BGFile)) {
          BCol <- c(BCol, grep(paste0("^", names(BGFile)[b], "$"), names(RandomOcc)))
        }
        RandomOcc <- RandomOcc[, c(BCol)]
        utils::write.csv(RandomOcc, file.path(model_output, spp.name, "NullModels", "RandomOcc.csv"),
                  row.names = FALSE)
        model.out <- tryCatch({
          #can turn off a lot of output writing for the final experiment
          #(jackknifing, write plot pngs) press help in maxent for details
          system(paste0("java -mx900m -jar ",  file.path(model_output, "maxent.jar"), " -e ", BackgroundFile, " -s ",
                        file.path(model_output, spp.name, "NullModels","RandomOcc.csv"),
                        " -J -o ", file.path(model_output, spp.name, "NullModels"),
                        " noaskoverwrite logistic threshold -X ",
                        test_percent, " replicates=", 1, " betamultiplier=", regularization,
                        " writeclampgrid=false", " writemess=false",
                        " nowarnings writeplotdata=false", " -a ",
                        reptype, " hinge=", hinge, " togglelayertype=", categorical))
        }, error = function(err) {
          print(paste("MY_ERROR: ", spp.name, " ", err))
          return(paste0("error: ", err))
        })
        failed_runs <- checkError(model.out, spp.name, failed_runs,
                                  paste0("model creation number ", i))
        if (!is.null(failed_runs)) {
          message(failed_runs)
        }
        MaxEntResults <- utils::read.csv(file.path(model_output, spp.name, "NullModels/maxentResults.csv"))
        NullAUCVals <- c(NullAUCVals, MaxEntResults$Test.AUC[1])
      }
      utils::write.csv(NullAUCVals, file.path(model_output, spp.name, "NullModel_AUC.csv"), row.names = FALSE)
      file.remove(model_output, "/", spp.name, "/NullModels")
    }
    
    for(i in 1:length(ListSpp)) {
      run(ListSpp[i])
    }
    
    # clus <- parallel::makeCluster(ncores, setup_timeout = 0.5)
    # 
    # envdata <- terra::wrap(envdata)
    # 
    # parallel::clusterExport(clus, varlist = c("run", "checkError", "occlist", "bglist", "ncores",
    #                                           "nrep", "categorical", "alloutputs", "model_output",
    #                                           "reptype","test_percent", "regularization", "hinge",
    #                                           "ListSpp", "replicates", "envdata", "bufflist"), envir = environment())
    # 
    # parallel::clusterEvalQ(clus, library(terra))
    # 
    # for (i in 1:nrow(ListSpp)) {
    #   out <- parallel::parLapply(clus, ListSpp[i, ], function(x) run(x))
    # }
    # 
    # parallel::stopCluster(clus)
  }

  do.call(NullMaxent, modelpar)

}

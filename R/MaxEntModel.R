#' Model species distributions with MaxEnt using parallel processing
#'
#' Takes occurrence points and background points of many species and models
#' them using the MaxEnt algorithm, parallelizing the process across multiple
#' computer cores.
#'
#' @param occlist a list of .csv file names, each containing the occurrence
#' points of a given taxon. The files should be named the same as the taxon
#' of interest (e.g.,: ".../Canis_lupus.csv").
#' @param bglist a list of .csv files corresponding to the background points
#' provided for each taxon. This list of files should be ordered in the same
#' way as the files given by \code{occlist}.
#' @param model_output the directory where all output files will be placed.
#' @param ncores the number of computer cores to parallelize the background point generation on.
#' Default is 1; Using one fewer core than the computer has is usually optimal.
#' @param nrep (integer) the number of replicates to run each species through.
#' @param categorical (character). If categorical variables are used for modelling (e.g., soil type),
#' they should be distinguished from the continuous data by a prefix (e.g., "C_soiltype.bil"). Provide
#' the distinguishing prefix here so that MaxEnt can distinguish bewteen categorical and continuous
#' environmental layers.
#' @param alloutputs Should secondary outputs from MaxEnt be generated \code{TRUE/FALSE}? Including:
#'
#'    1: A raster showing the spatial distribution of clamping for each run.
#'
#'    2: A multidimensional environmental similarity surface (MESS) showing novel climates.
#'
#'    3: Files containing the parameters used to make the response curves.
#'
#' The final set of arguments are optional and used for tuning the maxent model and cross-validation:
#' @param reptype Type of replication ("Crossvalidate", "Bootstrap", "Subsample"; see MAXENT manual).
#' Default is "Subsample".
#' @param test_percent (numeric): number between 0 and 100: percentage of points "held back" for
#' crossvalidation, Test AUC validation, etc. Default is 20.
#' @param features (optional): a vector of the features for MaxEnt to model the species-
#' environment relationships with. Options are one or more of \code{"linear", "quadratic", "product",
#' "threshold", "hinge"}. Refer to the MaxEnt help page for more information about each feature class
#' If there are few occurrence points, hinge features are discouraged. Default is all feature classes.
#' @param testsamples (optional) If cross-validation with a new set of occurrence points is required,
#' this should be a list of full file paths corresponding to the validation occurrence points for each
#' species. This will take presidence over the random test percentage given in \code{test_percent}.
#' NOTE: if using null AUC validation, testsamples must be given!
#' @param regularization (numeric) regularization parameter (penalizes complex models). A higher
#' regularization means more weight given to simpler models. Default is 1.
#' @export
#' @return Provides the trained model for each replicate and species (.lambdas file), a summary of
#' the outputs provided by the maxent.jar executable, a .csv file containing information on the
#' AUC values, threshold values, variable importance, etc., and (as requested) all of the outputs
#' given in the \code{alloutputs} description. A full summary of the output maxent.jar provides
#' can be found the MaxEnt manual.


MaxEntModel <- function(occlist, bglist, model_output,
                        ncores = 1, nrep = 1, categorical = NA,
                        alloutputs = TRUE, reptype = "Subsample",
                        test_percent = 20, features = c("linear", "quadratic", "product", "threshold", "hinge"),
                        testsamples = FALSE, regularization = 1) {

  if (length(occlist) < ncores) {
    ncores <- length(occlist)
  }
  ListSpp <- matrix(data = occlist, ncol = ncores)

  #Creates sub-directory for the given species
  
  if (!dir.exists(model_output)) {
    dir.create(model_output)
  }
  
  #Copies maxent.jar into the model_output folder
  if (!file.exists(file.path(model_output, "maxent.jar"))) {
    file.copy(from = system.file("extdata", "maxent.jar", package = "megaSDM"),
              to = file.path(model_output, "maxent.jar"))
  }
  
  if (!hasArg(features)) {
    linear <- "true"
    quadratic <- "true"
    product <- "true"
    threshold <- "true"
    hinge <- "true"
  } else {
    featuretypes <- c("linear", "quadratic", "product", "threshold", "hinge")
    for (i in 1:length(featuretypes)) {
      if (length(grep(featuretypes[i], features)) > 0) {
        assign(featuretypes[i], "true")
      } else {
        assign(featuretypes[i], "false")
      }
    }
  }

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

  run <- function(CurSpp) {
    s <- grep(CurSpp, occlist)

    OccurrenceFile <- occlist[s]
    BackgroundFile <- bglist[s]
    if (testsamples != FALSE) {
      TestFile <- testsamples[s]
    }
    #Determines the species name from the name of the occurrence file
    SpeciesSplit <- unlist(strsplit(OccurrenceFile, "/", fixed = TRUE))
    spp.name <- substr(SpeciesSplit[length(SpeciesSplit)], 1,
                          nchar(SpeciesSplit[length(SpeciesSplit)]) - 4)

    dir.create(file.path(model_output, spp.name))
    failed_runs <- c()

    if (is.logical(testsamples)) {
      #THIS IS THE MODEL CREATION COMMAND>>
      model.out <- tryCatch({
        #can turn off a lot of output writing for the final experiment
        #(jackknifing, write plot pngs) press help in maxent for details
        system(paste0("java -mx900m -jar ", file.path(model_output, "maxent.jar")," -e ", BackgroundFile, " -s ", OccurrenceFile,
                      " -J -o ", file.path(model_output, spp.name), " noaskoverwrite logistic threshold -X ",
                      test_percent, " replicates=", nrep, " betamultiplier=", regularization,
                      " writeclampgrid=", alloutputs, " writemess=", alloutputs,
                      " nowarnings writeplotdata=", alloutputs , " -a ",
                      reptype, " linear=", linear, " quadratic=", quadratic, " product=", product,
                      " threshold=", threshold, " hinge=", hinge, " togglelayertype=", categorical))
      }, error = function(err) {
        print(paste("MY_ERROR: ", spp.name, " ", err))
        return(paste0("error: ", err))
      })
    } else {
      model.out <- tryCatch({
        #can turn off a lot of output writing for the final experiment
        #(jackknifing, write plot pngs) press help in maxent for details
        system(paste0("java -mx900m -jar ", file.path(model_output, "maxent.jar")," -e ", BackgroundFile, " -s ", OccurrenceFile,
                      " -J -o ", file.path(model_output, spp.name), " noaskoverwrite logistic threshold -X ",
                      test_percent, " replicates=", 1, " betamultiplier=", regularization,
                      " testsamplesfile=", file.path(TestFile)," writeclampgrid=", alloutputs, " writemess=", alloutputs,
                      " nowarnings writeplotdata=", alloutputs , " -a ",
                      reptype, " linear=", linear, " quadratic=", quadratic, " product=", product,
                      " threshold=", threshold, " hinge=", hinge, " togglelayertype=", categorical))
      }, error = function(err) {
        print(paste("MY_ERROR: ", spp.name, " ", err))
        return(paste0("error: ", err))
      })
    }
    failed_runs <- checkError(model.out, spp.name, failed_runs,
                              paste0("model creation number ", i))
    if (!is.null(failed_runs)) {
      message(failed_runs)
    }

  }

  clus <- parallel::makeCluster(ncores, setup_timeout = 0.5)

  parallel::clusterExport(clus, varlist = c("run", "checkError", "occlist", "bglist", "ncores",
                                            "nrep", "categorical", "alloutputs", "model_output",
                                            "reptype","test_percent", "regularization", "linear",
                                            "quadratic", "product", "threshold", "hinge",
                                            "testsamples", "ListSpp"), envir = environment())

  parallel::clusterEvalQ(clus, library(raster))

  for (i in 1:nrow(ListSpp)) {
    out <- parallel::parLapply(clus, ListSpp[i, ], function(x) run(x))
  }

  parallel::stopCluster(clus)
}

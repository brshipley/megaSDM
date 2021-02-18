context("Environmental Data Management")
library(megaSDM)
library(raster)
test_that("All example environmental data can be loaded", {
  Train <- list.files(system.file("extdata", "trainingarea", package = "megaSDM"),
                      pattern = ".bil$",
                      full.names = TRUE)
  Proj <- list.files(system.file("extdata", "predictenv/RCP4.5/2050", package = "megaSDM"),
                     pattern = ".bil$",
                     full.names = TRUE)
  expect_output(raster(Train[1]), "RasterLayer")
  expect_outpu(raster(Proj[1]), "RasterLayer")
})

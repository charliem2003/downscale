################################################################################
# 
# Hui.downscale.R
# Version 1.3
# 27/07/2017
#
# Updates:
#   27/07/2017: SpatialPointsDataFrame allowed as input
#               'lat' and 'lon' as column names replaced by 'x' and 'y'
#   29/09/2016: No need to input cell width in hui.downscale with raster
#   05/05/2014: calculates AOO
#               Can use upgrain as an input
# 
# Functions to calculate q++ - the conditional probability that a randomly
# chosen cell adjacent to an occupied cell is also occupied
#
# Args:
#   atlas.data: either a data frame containing a column of presence and absence 
#               data and a column each for longitude and latitude; or a raster
#               file where 1 = presence and 0 = absence; or a 
#               SpatialPointsDataFrame with a data frame containing 'presence'
#   cell.width: Cell area of the atlas data (ie resolution)
#   new.areas: vector of grain sizes for model prediction (area)
#   extent: extent of points for conversion to AOO
#   tolerance: tolerance for root solving to estimate probability of absence at
#              the fine scale
#   plot: pass to plot.predict.downscale
#
################################################################################

hui.downscale <- function(atlas.data,
                          cell.width,
                          new.areas,
                          extent = NULL,
                          tolerance = 1e-6,
                          plot = FALSE) {
  
  ### error checking: if data frame requires extent
  if(is.data.frame(atlas.data) == "TRUE") {
    if(is.null(extent)) {
      stop("Extent required if data input is data frame of coordinates")
    }
  }
  
  ### error checking: if SpatialPointsDataFrame requires extent
  if(class(atlas.data)[1] == "SpatialPointsDataFrame") {
    if(is.null(extent)) {
      stop("Extent required if data input is SpatialPointsDataFrame")
    }
  }
  
  ### Error checking: if data frame needs cell width
  if(is.data.frame(atlas.data) == "TRUE") {
    if(is.null(cell.width)) {
      stop("If data is data.frame cell.width is required")
    }
  }
  
  ### Error checking: if SpatialPointsDataFrame needs cell width
  if(class(atlas.data)[1] == "SpatialPointsDataFrame") {
    if(is.null(cell.width)) {
      stop("If data is SpatialPointsDataFrame cell.width is required")
    }
  }
  
  ###############################################################################
  ### data input handling
  if(class(atlas.data) == "upgrain") {
    extent <- atlas.data$extent.stand
    species <- raster::rasterToPoints(atlas.data$atlas.raster.stand)
    species <- data.frame(presence = species[, 3],
                          x = species[, "y"],
                          y = species[, "x"])
    cell.width <- raster::res(atlas.data$atlas.raster.stand)[1]
    cell.area <- cell.width ^ 2
  }
  
  if(class(atlas.data) == "RasterLayer") {
    species <- raster::rasterToPoints(atlas.data)
    species <- data.frame(presence = species[, 3],
                          x = species[, "y"],
                          y = species[, "x"])
    cell.width <- raster::res(atlas.data)[1]
    cell.area <- cell.width ^ 2
  }
  
  if(class(atlas.data) == "SpatialPointsDataFrame") {
    species <- data.frame(presence = atlas.data@data[, "presence"],
                          x = atlas.data@coords[, "y"],
                          y = atlas.data@coords[, "x"])
    cell.area <- cell.width ^ 2
  }
  
  if(is.data.frame(atlas.data) == "TRUE") {
    species <- atlas.data
    cell.area <- cell.width ^ 2
  }
  
  ##############################################################################
  ### error checking that it's presence-absence data
  if(sum(species[, "presence"] > 1, na.rm = TRUE) > 0) {
    stop("Presence-absence data not in 1's and 0's")
  }
  
  ### error checking: predicting at finer scales than atlas scale
  if(sum(new.areas >= cell.area) > 0) {
    stop("One or more fine scale grid sizes are larger than atlas scale grid size",
         call. = FALSE)
  }
  
  #############################################################################
  ### Hui modelling
  
  # p+ : observed probability of occurrence at coarse grain
  p1_coarse <- sum(species[, "presence"] == 1, na.rm = TRUE) / 
    sum(!is.na(species[, "presence"]))
  
  # p- : observed probability of absence at coarse grain
  p0_coarse <- 1 - p1_coarse  
  
  # q+/+ :  conditional probability that a randomly chosen cell adjacent to an 
  #         occupied cell is also occupied
  q11 <- prob_q11(presences = species, unit = cell.width)
  
  # the conditional probability that a randomly chosen cell adjacent to an 
  # empty cell is occupied e.g. 1 - q+/0 at the coarse grain
  q00 <- (1 - (2 * p1_coarse) + (q11 * p1_coarse)) / (1 - p1_coarse)
  
  # solve for p0_fine - probability of absence at fine scale
  predicted <- data.frame("Cell.area" = new.areas, "Occupancy" = NA, "AOO" = NA)
  for(i in 1:length(new.areas)) {
    fine.width <- sqrt(new.areas[i])
    p0_fine_root <- uniroot(f = ResidHui,
                       n = cell.width / fine.width,
                       q00 = q00,
                       p0_coarse = p0_coarse,
                       lower= 1e-12,
                       upper = 1,
                       tol = tolerance,
                       extendInt = "yes")
    predicted[i, "Occupancy"] <- 1 - p0_fine_root$root
  }
  predicted[, "AOO"] <- predicted[, "Occupancy"] * extent
  observed <- data.frame("Cell.area" = cell.width ^ 2,
                         "Occupancy" = p1_coarse)
  output <- list("model" = "Hui",
                 "predicted" = predicted,
                 "observed" = observed)
  class(output) <- "predict.downscale"
  
  if (plot == TRUE) {
    par.original <- par()
    par.original <- list(mfrow = par.original$mfrow, mar = par.original$mar)
    par(mfrow = c(1, 1), mar = c(5, 5, 3, 1))
    
    plot.predict.downscale(output)
    
    par(mfrow = par.original$mfrow, mar = par.original$mar)
  }
  return(output)
}
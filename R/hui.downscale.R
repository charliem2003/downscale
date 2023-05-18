################################################################################
# 
# Hui.downscale.R
# Version 2.0
# 19/05/2023
#
# Updates:
#   19/05/2023: v2.0 - CONVERTED TO TERRA AND SF
#               error checking moved to checkInputs
#   26/10/2021: uses on.exit to return to original par settings
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
#               file where 1 = presence and 0 = absence; or an sf spatial points
#               object with a data frame containing 'presence'.
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
  
  ##############################################################################
  ### Error checking
  checkInputs(inputFunction = "hui.downscale",
              atlas.data = atlas.data,
              cell.width = cell.width,
              extent = extent) 
  
  ###############################################################################
  ### data input handling
  if(class(atlas.data)[1] %in% c("upgrain", "SpatRaster")) {
    cell.area <- cell.width ^ 2
    
    if(class(atlas.data)[1] == "upgrain") {
      species <- atlas.data$atlas.raster.stand
      names(species)[1] <- "presence"
      cell.width <- terra::res(species)[1]
      extent  <- atlas.data$extent.stand
    }
    
    if(class(atlas.data)[1] == "SpatRaster") {
      species <- atlas.data
      names(species)[1] <- "presence"
      cell.width <- terra::res(species)[1]
      extent  <- sum(!is.na(terra::values(species))) * (cell.width ^ 2)
    }
    
    cells <- which(!is.na(values(species)))
    species <- data.frame(presence = species[cells],
                          x        = terra::xFromCell(species, cells),
                          y        = terra::yFromCell(species, cells))
  }
  
  if(class(atlas.data)[1] == "sf") {
    species <- data.frame(presence = atlas.data$presence,
                          x =        sf::st_coordinates(atlas.data)[, "X"],
                          y =        sf::st_coordinates(atlas.data)[, "Y"])
  }

  if(class(atlas.data)[1] == "data.frame") {
    species <- atlas.data
  }
  cell.area <- cell.width ^ 2
  
  ##############################################################################
  ### Error checking: that it's presence-absence data
  if(sum(species$presence > 1, na.rm = TRUE) > 0) {
    stop("Presence-absence data not in 1's and 0's", call. = FALSE)
  }
  
  ### Error checking: predicting at finer scales than atlas scale
  if(sum(new.areas >= cell.area) > 0) {
    stop("One or more fine scale grid sizes are larger than atlas scale grid size",
         call. = FALSE)
  }
  
  #############################################################################
  ### Hui modelling
  
  # p+ : observed probability of occurrence at coarse grain
  p1_coarse <- sum(species$presence == 1, na.rm = TRUE) /
    sum(!is.na(species$presence))
  
  # p- : observed probability of absence at coarse grain
  p0_coarse <- 1 - p1_coarse  
  
  # q+/+ :  conditional probability that a randomly chosen cell adjacent to an 
  #         occupied cell is also occupied
  q11 <- prob_q11(presences = species, unit = cell.width)
  
  # the conditional probability that a randomly chosen cell adjacent to an 
  # empty cell is occupied e.g. 1 - q+/0 at the coarse grain
  q00 <- (1 - (2 * p1_coarse) + (q11 * p1_coarse)) / (1 - p1_coarse)
  
  # solve for p0_fine - probability of absence at fine scale
  predicted <- data.frame(Cell.area = new.areas,
                          Occupancy = NA,
                          AOO       = NA)
  
  for(i in 1:length(new.areas)) {
    fine.width <- sqrt(new.areas[i])
    p0_fine_root <- uniroot(f         = ResidHui,
                            n         = cell.width / fine.width,
                            q00       = q00,
                            p0_coarse = p0_coarse,
                            lower     = 1e-12,
                            upper     = 1,
                            tol       = tolerance,
                            extendInt = "yes")
    predicted$Occupancy[i] <- 1 - p0_fine_root$root
  }
  
  ### output object
  predicted$AOO <- predicted$Occupancy * extent
  observed <- data.frame(Cell.area = cell.width ^ 2,
                         Occupancy = p1_coarse)
  output <- list(model     = "Hui",
                 predicted = predicted,
                 observed  = observed)
  class(output) <- "predict.downscale"
  
  ### optional plotting
  if (plot == TRUE) {
    parOrig <- par(no.readonly = TRUE)
    on.exit(par(parOrig))
    par(mfrow = c(1, 1), mar = c(5, 5, 3, 1))
    
    plot.predict.downscale(output)
  }
  return(output)
}

################################################################################
# 
# upgrain.R
# Version 2.0
# 19/05/2023
#
# Updates:
#   19/05/2023: v2.0 - CONVERTED TO TERRA AND SF
#               error checking moved to checkInputs
#   16/05/2023: Simple reformatting
#   26/10/2021: uses on.exit to return to original par settings
#   31/07/2017: Bug when using SpatialPointsDataFrame fixed
#   27/07/2017: SpatialPointsDataFrame allowed as input
#               'lat' and 'lon' as column names replaced by 'x' and 'y'
#   18/08/2016: Change default to "All_Sampled"
#   25/02/2016: Bug on extent in standardised rasters fixed
#               Added option to save all rasters
#               Change the method of assigning raster values to 'setValues'
#   05/10/2015: Bug on number of scales fixed
#               'All_Presences' changed to 'All_occurrences'
#               Plotting now optional
#   30/04/2015: Threshold selection added
#   14/04/2015: plotting functions added
#               help file updated
#
# Takes presence-absence atlas data and aggregates data to larger grain sizes, 
# returning occupancy at each grain size for use in downscale modelling. The 
# extent for all scales is standardised to that of the largest grain size by 
# applying a threshold for the proportion of unsampled atlas cells allowed 
# within a cell at the largest grain size. The threshold can be chosen by 
# the user or one of four threshold selections methods apply. 
#
# # args:
#   atlas.data = either a raster file of presence-absence atlas data, or a data 
#                frame of sampled cells with longitude, latitude and 
#                presence-absence; or a sf spatial points object with a data 
#                frame containing 'presence'.
#   cell.width = if data is a data frame, the cell widths of sampled cells.
#   scales = the number of cells to upgrain. Upgraining will happen by factors 
#            of 2 - ie if scales = 3, the atlas data will be aggregated in 2x2 
#            cells, 4x4 cells and 8x8 cells.
#   threshold = a value for the proportion of unsampled atlas cells allowed 
#               within a cell at the largest grain size. Default = NULL.
#   method = one of the default methods for selecting a threshold. There are
#            four choices: "All_Sampled" (default), "All_Occurrences", 
#            "Gain_Equals_Loss" or "Sampled_Only".
#   plot = Plots the original atlas data alongside the standardised atlas data 
#          at each grain size. Default = TRUE.
#   return.rasters = If TRUE returns the extent-standardised atlas data 
#                    upgrained to all grain sizes (NOTE: the extent-standardised
#                    atlas data at the original grain size is always returned
#                    regardless). Default = FALSE.
#
# outputs:
#   An object of class 'upgrain', containing five objects:
#     threshold = a value for the threshold used.
#     extent.stand = a value for the standardised extent.
#     occupancy.stand = data frame of cell area, extent and occupancies of after
#                       upgrain-standardisation to be used as input in
#                       downscaling
#     occupancy.orig = data frame of cell area, extent and occupancies before
#                      upgrain-standardisation.
#     atlas.raster.stand = A raster layer of the extent-standardised atlas data.
#     scaled.rasters = (if return.rasters = TRUE) a list containing the extent-
#                      standardised atlas data upgrained to all grain sizes
#
################################################################################

upgrain <- function(atlas.data,
                    cell.width = NULL,
                    scales,
                    threshold = NULL,
                    method = "All_Sampled",
                    plot = TRUE,
                    return.rasters = FALSE) {
  
  ##############################################################################
  ### Error checking
  checkInputs(inputFunction = "upgrain",
                       atlas.data    = atlas.data,
                       cell.width    = cell.width,
                       scales        = scales,
                       threshold     = threshold,
                       method        = method)
  
  ##############################################################################
  ### data storage
  original <- data.frame(Cell.area = rep(NA, scales + 1), 
                         Extent    = rep(NA, scales + 1),
                         Occupancy = rep(NA, scales + 1))
  extended <- data.frame(Cell.area = rep(NA, scales + 1), 
                         Extent    = rep(NA, scales + 1),
                         Occupancy = rep(NA, scales + 1))
  
  ##############################################################################
  ### data manipulation - data.frame -> sf -> raster (depending on input data)
  
  ### if data frame convert to spatial points first
  if(class(atlas.data)[1] == "data.frame") {
    atlas.data$presence[atlas.data$presence > 0] <- 1
    atlas.data <- sf::st_as_sf(atlas.data, coords = c("x", "y"))
  }
  
  ### if sp object convert to sf
  if(inherits(atlas.data, "SpatialPointsDataFrame")) {
    warning("The 'sp' package is becoming obsolete.
In future downscale updates 'SpatialPointsDataFrame' objects will be phased out.
We recommend using the 'sf' package instead",
         call. = FALSE)
    atlas.data <- sf::st_as_sf(atlas.data)
  }
  
  ### rasterize spatial points
  if(class(atlas.data)[1] == "sf") {
    ext <- sf::st_bbox(atlas.data)
    atlas_raster <- terra::rast(xmin = ext$xmin - (cell.width / 2),
                                xmax = ext$xmax + (cell.width / 2),
                                ymin = ext$ymin - (cell.width / 2),
                                ymax = ext$ymax + (cell.width / 2),
                                resolution = cell.width)
    atlas_raster$presence <- terra::rasterize(atlas.data, atlas_raster,
                                              field = "presence", fun = max)
  }
  
  ### if raster object convert to terra
  if(inherits(atlas.data, "RasterLayer")) {
    warning("The 'raster' package is becoming obsolete.
In future downscale updates using 'raster' class objects will be phased out.
We recommend using a 'SpatRaster' object using the 'terra' package instead",
            call. = FALSE)
    atlas.data <- terra::rast(atlas.data)
  }
  
  ### if raster object
  if(class(atlas.data)[1] == "SpatRaster") {
    atlas_raster <- atlas.data
  }
  
  ### check that only 1s and 0s
  atlas_raster[atlas_raster > 0] <- 1
  
  ### largest scale (all other layers are extended to equal this raster)
  max_raster <- terra::aggregate(atlas_raster, fact = (2 ^ scales),
                                 fun = "max", na.rm = TRUE)
  atlas_raster_extend <- terra::extend(atlas_raster, terra::ext(max_raster))
  
  ##############################################################################
  ### Set new atlas data to selected threshold (atlas_thresh)
  
  ### create boundary raster
  boundary_raster <- atlas_raster
  boundary_raster[!is.na(boundary_raster)] <- 1
  boundary_raster <- terra::aggregate(boundary_raster, fact = (2 ^ scales),
                                      fun = "sum", na.rm = TRUE)
  boundary_raster <- boundary_raster / ((2 ^ scales) ^ 2)
  
  ### if no threshold specified determine threshold base on method
  if(is.null(threshold)) {
    ### set thresholds for "Sampled_Only" and "All_Sampled"
    if(method == "Sampled_Only") {
      threshold <- 1
    }
    if(method == "All_Sampled") {
      threshold <- 0
    }
  
    ### determine thresholds for "All_Occurrences" or "Gain_Equals_Loss" methods
    if(method %in% c("All_Occurrences", "Gain_Equals_Loss")) {
      atlas_boundary <- atlas_raster_extend
      atlas_boundary[atlas_boundary == 0]   <- 1
      atlas_boundary[is.na(atlas_boundary)] <- 0
      
      ### Loop through potential thresholds
      thresholds <- seq(0, 1, 0.01)
      
      land <- data.frame(Threshold           = thresholds, 
                         SampledExcluded     = rep(NA, length(thresholds)),
                         SampledIncluded     = NA,
                         UnsampledAdded      = NA, 
                         Extent              = NA,
                         OccurrencesExcluded = NA)
      
      for(j in 1:length(thresholds)) {
        ### redraw boundary based on threshold
        atlas_boundary_thresh <- max_raster
        atlas_boundary_thresh[boundary_raster < thresholds[j]] <- NA
        atlas_boundary_thresh[atlas_boundary_thresh >= 0]   <- 2
        atlas_boundary_thresh[is.na(atlas_boundary_thresh)] <- 0
        atlas_boundary_thresh <- terra::disagg(atlas_boundary_thresh,
                                               fact = (2 ^ scales),
                                               method = "near")
        both <- sum(atlas_boundary_thresh, atlas_boundary)
        
        ### store values in data frame
        land$SampledExcluded[j] <- sum(terra::values(both) == 1, na.rm = TRUE)
        land$UnsampledAdded[j]  <- sum(terra::values(both) == 2, na.rm = TRUE)
        land$SampledIncluded[j] <- sum(terra::values(both) == 3, na.rm = TRUE)
        land$Extent[j]          <- sum(terra::values(atlas_boundary_thresh) ==2,
                                       na.rm = TRUE)
        
        ### proportion occurrences that fall outside boundary
        both <- sum(atlas_boundary_thresh, atlas_raster_extend)
        excl <- sum(terra::values(both) == 1, na.rm = TRUE) /
                sum(terra::values(atlas_raster) == 1, na.rm = TRUE)
        land$OccurrencesExcluded[j] <- round(excl, 3)
      }
      
      ### assign appropriate thresholds
      if(method == "All_Occurrences") {
        threshold <- thresholds[max(which(land[, "OccurrencesExcluded"] == 0))]
      }
      if(method == "Gain_Equals_Loss") {
        diffs <- land$Extent - sum(terra::values(atlas_boundary), na.rm = TRUE)
        threshold <- thresholds[which.min(abs(diffs))]
      }
    }
  }
  
  ### if threshold is still NULL then there is a problem
  if(is.null(threshold)) {
    stop("Unable to determine a threshold and no threshold suggested")
  }
  
  ### set extent of maximum grain size raster to threshold (max_raster_thresh)
  max_raster_thresh <- max_raster
  max_raster_thresh[boundary_raster < threshold] <- NA
  max_raster_thresh[!is.na(max_raster_thresh)] <- 1
  max_raster_thresh <- terra::disagg(max_raster_thresh,
                                     fact = (2 ^ scales), method = "near")
  
  ### extend and crop atlas raster to new extent (atlas_thresh)
  atlas_thresh <- atlas_raster
  atlas_thresh <- terra::extend(atlas_thresh, terra::ext(max_raster_thresh))
  atlas_thresh[atlas_thresh == 1]   <- 2
  atlas_thresh[atlas_thresh == 0]   <- 1
  atlas_thresh[is.na(atlas_thresh)] <- 0
  
  atlas_thresh <- sum(atlas_thresh, max_raster_thresh)
  atlas_thresh[atlas_thresh == 1] <- 0
  atlas_thresh[atlas_thresh == 2] <- 0
  atlas_thresh[atlas_thresh == 3] <- 1
  
  ##############################################################################
  #### calculate occupancy at all grain sizes in original atlas data
  
  original$Cell.area[1] <- terra::res(atlas_raster)[1] ^ 2
  original$Extent[1]    <- sum(!is.na(terra::values(atlas_raster))) *
                           original$Cell.area[1]
  original$Occupancy[1] <- sum(terra::values(atlas_raster) == 1, na.rm = TRUE) /
                           sum(!is.na(terra::values(atlas_raster)))
  
  for(i in 1:scales) {
    scaled_raster <- terra::aggregate(atlas_raster, fact = (2 ^ i), fun = "max")
    original$Cell.area[i + 1] <- terra::res(scaled_raster)[1] ^ 2
    original$Extent[i + 1]    <- sum(!is.na(terra::values(scaled_raster))) *
                                 original$Cell.area[i + 1]
    original$Occupancy[i + 1] <- sum(terra::values(scaled_raster) == 1,
                                   na.rm = TRUE) /
                                 sum(!is.na(terra::values(scaled_raster)))
  }
  
  ##############################################################################
  #### calculate occupancy at all grain sizes in standardised atlas data
  
  extended$Cell.area[1] <- terra::res(atlas_thresh)[1] ^ 2
  extended$Extent[1]    <- sum(!is.na(terra::values(atlas_thresh))) * 
                           extended$Cell.area[1]
  extended$Occupancy[1] <- sum(terra::values(atlas_thresh == 1, na.rm = TRUE)) /
                           sum(!is.na(terra::values(atlas_thresh)))
  
  for(i in 1:scales) {
    scaled_raster <- terra::aggregate(atlas_thresh, fact = (2 ^ i), fun = "max")
    if(return.rasters == TRUE) {
      assign(paste("scaled_raster", i, sep = "_"), scaled_raster)
    }
    extended$Cell.area[i + 1] <- terra::res(scaled_raster)[1] ^ 2
    extended$Extent[i + 1]    <- sum(!is.na(terra::values(scaled_raster))) *
                                 extended$Cell.area[i + 1]
    extended$Occupancy[i + 1] <- sum(terra::values(scaled_raster) == 1,
                                        na.rm = TRUE) /
                                 sum(!is.na(terra::values(scaled_raster)))
  }
  extended$Cell.area <- original$Cell.area
  
  ##############################################################################
  ### plotting
  
  if(plot == TRUE) {  
    ### set par values for plotting
    parOrig <- par(no.readonly = TRUE)
    on.exit(par(parOrig))
    par(mfrow = c(2, ceiling((scales + 2) / 2)), mar = c(3, 3, 3.5, 3))
    
    ########### Plot 1 - original atlas data
    terra::plot(atlas_raster_extend,
                col = c("white", "white"),
                legend = FALSE, 
                axes = FALSE,
                mar = c(3, 3, 3.5, 3),
                main = paste0("Original atlas data:\n cell area = ", 
                              original$Cell.area[1]))
    terra::plot(atlas_raster,
                col = c("white", "red"),
                legend = FALSE, 
                axes = FALSE,
                mar = c(3, 3, 3.5, 3),
                colNA = "dark grey",
                add = TRUE)
    rect(xleft   = terra::ext(atlas_raster)$xmin,
         ybottom = terra::ext(atlas_raster)$ymin,
         xright  = terra::ext(atlas_raster)$xmax,
         ytop    = terra::ext(atlas_raster)$ymax)
    
    rect(xleft   = terra::ext(atlas_raster_extend)$xmin,
         ybottom = terra::ext(atlas_raster_extend)$ymin,
         xright  = terra::ext(atlas_raster_extend)$xmax,
         ytop    = terra::ext(atlas_raster_extend)$ymax, lty = 2)
    
    ########### Plot 2 - standardised atlas data - atlas grain size
    terra::plot(atlas_thresh,
                col = c("white", "red"),
                legend = FALSE, 
                axes = FALSE,
                mar = c(3, 3, 3.5, 3),
                colNA = "dark grey",
                main = paste0("Standardised atlas data:\n cell area = ", 
                              original$Cell.area[1]))
    rect(xleft   = terra::ext(atlas_thresh)$xmin,
         ybottom = terra::ext(atlas_thresh)$ymin,
         xright  = terra::ext(atlas_thresh)$xmax,
         ytop    = terra::ext(atlas_thresh)$ymax)
    
    for(i in 1:scales) {
      scaled_raster <- terra::aggregate(atlas_thresh,
                                        fact = (2 ^ i),
                                        fun = "max")
      
      ########### Plots 3 ... - standardised atlas data - all other grain sizes
      terra::plot(scaled_raster,
                  col = c("white", "red"),
                  colNA = "dark grey",
                  legend = FALSE, 
                  axes = FALSE,
                  mar = c(3, 3, 3.5, 3),
                  main = paste0("Standardised atlas data:\n cell area = ", 
                                original$Cell.area[i + 1]))
      rect(xleft   = terra::ext(scaled_raster)$xmin,
           ybottom = terra::ext(scaled_raster)$ymin,
           xright  = terra::ext(scaled_raster)$xmax,
           ytop    = terra::ext(scaled_raster)$ymax)
    }
  }
  
  ##############################################################################
  ### output
  
  output <- list(threshold          = threshold,
                 extent.stand       = extended$Extent[1],
                 occupancy.stand    = extended,
                 occupancy.orig     = original,
                 atlas.raster.stand = atlas_thresh)
  
  if(return.rasters == TRUE) {
    scaled.rasters <- list()
    for(i in 1:scales) {
      scaled.rasters[i] <- get(paste("scaled_raster", i, sep = "_"))
    }
    names(scaled.rasters) <- paste0("scaled.raster", extended$Cell.area[-1])
    output <- list(threshold          = threshold,
                   extent.stand       = extended[1, "Extent"],
                   occupancy.stand    = extended,
                   occupancy.orig     = original,
                   atlas.raster.stand = atlas_thresh,
                   scaled.rasters     = scaled.rasters)
  }
  
  class(output) <- "upgrain"
  return(output)
}
################################################################################
#
# upgrain.threshold.R 
# Version 2.0
# 19/05/2023
#
# Updates:
#   19/05/2023: v2.0 - CONVERTED TO TERRA AND SF
#               slight changes to plots
#               error checking moved to checkInputs
#   16/05/2023: Simple reformatting
#   15/11/2021: plotting option added;
#               bugfix when boundary_poly is null
#   26/10/2021: uses on.exit to return to original par settings
#   31/07/2017: Bug when using SpatialPointsDataFrame fixed
#   27/07/2017: SpatialPointsDataFrame allowed as input
#               'lat' and 'lon' as column names replaced by 'x' and 'y'
#   05/10/2015: Bug on number of scales fixed
#               'All_Presences' changed to 'All_occurrences'
#
# Explores the NoData threshold selection of cell selection for upgraining. A 
# threshold of 0 means that all cells at the largest grain with some Sampled 
# are included (All Sampled). A threshold of 1 means only cells at the largest 
# grain that are entirely within the Sampled are included (Sampled only). Two 
# other threshold options are also suggested based upon the trade-off of 
# including NoData cells as absence and losing presence records: 
# All Occurrences: the threshold that preserves all presence records
# Gain Equals Loss: the threshold where the upgrained extent = the original 
# extent
#
# args:
#   atlas.data = either a raster file of presence-absence atlas data, or a data 
#                frame of sampled cells with longitude, latitude and 
#                presence-absence; or an sf spatial points object with a data 
#                frame containing 'presence'.
#   cell.width = if data is a data frame, the cell widths of sampled cells.
#   scales     = the number of cells to upgrain. Upgraining will happen by
#                factors of 2 - ie if scales = 3, the atlas data will be
#                aggregated in 2 x 2 cells, 4 x 4 cells and 8 x 8 cells.
#   thresholds = a vector of thresholds between and including 0 and 1 for the
#                quantity of Unsampled NoData cells that can be included.
#   plot       = See outputs section for plot details. Default = TRUE.
#
# outputs:
#   Two plotting windows. One with four plots relevant to exploring the 
#   trade-offs. The second plots the maps of the four default threshold 
#   selections.
#   A list with two objects;
#     Thresholds: the threshold values for the four default threshold selections
#     Data: Data frame with all the relevant data for exploring the trade-offs.
#
################################################################################

upgrain.threshold <- function(atlas.data,
                              cell.width  = NULL,
                              scales,
                              thresholds = seq(0, 1, 0.01),
                              plot = TRUE) {
 
  ##############################################################################
  ### Error checking
  checkInputs(inputFunction = "upgrain.threshold",
              atlas.data    = atlas.data,
              cell.width    = cell.width,
              scales        = scales,
              thresholds    = thresholds)
  
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
  
  ########################################################
  ### create boundary raster
  boundary_raster <- atlas_raster
  boundary_raster[!is.na(boundary_raster)] <- 1
  boundary_raster <- terra::aggregate(boundary_raster, fact = (2 ^ scales),
                                      fun = "sum", na.rm = TRUE)
  boundary_raster <- boundary_raster / ((2 ^ scales) ^ 2)
  boundary_poly <- terra::as.polygons(boundary_raster, dissolve = TRUE)
  
  ##################################################
  ### Loop through thresholds
  thresholds <- thresholds
  land <- data.frame(Threshold           = thresholds, 
                     SampledExcluded     = rep(NA, length(thresholds)),
                     SampledIncluded     = NA,
                     UnsampledAdded      = NA, 
                     Extent              = NA,
                     OccurrencesExcluded = NA)
  
  atlas_raster_extend <- terra::extend(atlas_raster, terra::ext(max_raster))
  atlas_boundary <- atlas_raster_extend
  atlas_boundary[atlas_boundary == 0]   <- 1
  atlas_boundary[is.na(atlas_boundary)] <- 0
  
  for(j in 1:length(thresholds)){
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
  
  ### Possible Thresholds
  diffs <- land$Extent - sum(terra::values(atlas_boundary), na.rm = TRUE)
  gain_loss_thresh <- thresholds[which.min(abs(diffs))]
  
  occ_thresh <- thresholds[max(which(land[, "OccurrencesExcluded"] == 0))]
  
  ##############################################################################
  ### Plotting
  
  if(plot == TRUE) {  
    parOrig <- par(no.readonly = TRUE)
    on.exit(par(parOrig))
    par(mfrow = c(2, 2), mar = c(5, 5, 1, 1))
    
    ########################################
    ### First set of plots - diagnostics ###
    
    ### Plot 1: the total extent
    plot(land$Threshold, land$Extent, type = "l",
         ylab = "Extent (total number of cells)", xlab = "Threshold",
         ylim = c(0, max(land$Extent, na.rm = TRUE)))
    abline(h = sum(values(atlas_boundary) == 1, na.rm = TRUE), col = "Grey")
    lines(c(gain_loss_thresh, gain_loss_thresh),
          c(0, land$Extent[thresholds == gain_loss_thresh]),
          col = "red", lty = 2)
    lines(c(occ_thresh, occ_thresh),
          c(0, land$Extent[thresholds == occ_thresh]),
          col = "blue", lty = 2)
    legend("topright",
           legend = c("Original extent",
                      "Gain equals loss threshold",
                      "All occurrences threshold"),
           lty = c(1, 2, 2), bty = "n", col = c("Grey", "Red", "Blue"))
    
    ### Plot 2: the number of Sampled cells incorrectly identified
    plot(land$Threshold, land$SampledExcluded, type = "l",
         ylab = "Number of cells", xlab = "Threshold",
         ylim = c(0, max(c(land$SampledExcluded, land$UnsampledAdded),
                         na.rm = TRUE)))
    lines(land[, "Threshold"], land[, "UnsampledAdded"], lty = 2)
    lines(c(gain_loss_thresh, gain_loss_thresh),
          c(0, land$SampledExcluded[thresholds == gain_loss_thresh]),
          col = "red", lty = 2)
    lines(c(occ_thresh, occ_thresh),
          c(0, land$UnsampledAdded[thresholds == occ_thresh]),
          col = "blue", lty = 2)
    legend("top",
           legend = c("Number of sampled cells excluded",
                      "Number of unsampled cells added"),
           lty = 1:2, bty = "n")
    
    ### Plot 3: the number of original land cells retained
    plot(land$Threshold,
         land$SampledIncluded / max(land$SampledIncluded, na.rm = TRUE),
         type = "l",
         ylab = "Prop. of sampled cells retained", xlab = "Threshold",
         ylim = c(0, 1))
    lines(c(gain_loss_thresh, gain_loss_thresh),
          c(0, (land$SampledIncluded[thresholds == gain_loss_thresh] /
                  max(land$SampledIncluded, na.rm = TRUE))),
          col = "red", lty = 2)
    lines(c(occ_thresh, occ_thresh),
          c(0, (land$SampledIncluded[thresholds == occ_thresh] /
                  max(land$SampledIncluded, na.rm = TRUE))),
          col = "blue", lty = 2)
    
    ### Plot 4: % of presence cells excluded
    plot(land$Threshold, land$OccurrencesExcluded, type = "l",
         ylab = "Prop. of occurrences excluded", xlab = "Threshold",
         ylim = c(0, 1))
    lines(c(gain_loss_thresh, gain_loss_thresh),
          c(1, land$OccurrencesExcluded[thresholds == gain_loss_thresh]),
          lty = 2, col = "red")
    lines(c(occ_thresh, occ_thresh),
          c(1, land$OccurrencesExcluded[thresholds == occ_thresh]),
          lty = 2, col = "blue")
    
    ########################################################
    ### Second set of plots - Maps of thresholds options ###
    
    par(ask = TRUE) 
    par(mfrow = c(2, 2), mar = c(5.5, 1, 3.5, 1))
    
    thresh.selection <- c(0, occ_thresh, gain_loss_thresh, 1)
    selection <- c("All Sampled",
                   "All Occurrences",
                   "Gain Equals Loss",
                   "Sampled Only")
    selection <- selection[order(thresh.selection)]
    thresholds <- thresholds[order(thresh.selection)] 
    
    for(j in 1:length(thresh.selection)) {
      thresh <- thresh.selection[j]
      
      ### raster of selected cells
      max_raster_thresh <- max_raster
      max_raster_thresh[boundary_raster < thresh] <- NA
      
      ### polygon of thresholded raster boundary
      max_raster_thresh_extent <- max_raster_thresh
      max_raster_thresh_extent[max_raster_thresh_extent == 0] <- 1
      boundary_poly <- terra::as.polygons(max_raster_thresh_extent,
                                          dissolve = TRUE)
      
      # max_raster_thresh <- terra::disagg(max_raster_thresh, fact = 8)
      
      # ### extend original atlas data to extent of thresholded raster
      # atlas_boundary <- atlas_raster
      # atlas_boundary[atlas_boundary == 0]   <- 1
      # atlas_boundary[is.na(atlas_boundary)] <- 0
      # atlas_boundary <- terra::extend(atlas_boundary,
      #                                 terra::ext(max_raster_thresh))
      
      ### final thresholded atlas
      atlas_boundary_thresh <- max_raster_thresh
      atlas_boundary_thresh[atlas_boundary_thresh >= 0] <- 2
      atlas_boundary_thresh[is.na(atlas_boundary_thresh)] <- 0
      
      ### map of original atlas data extended to full extent
      terra::plot(atlas_raster_extend, axes = FALSE,
                  colNA = rgb(0.5, 0.5, 0.5),
                  col = c(rgb(1, 1, 1), rgb(1, 0, 0)),
                  legend = FALSE,
                  mar = c(5.5, 1, 3.5, 1))
      title(main = paste(selection[j], "\n", "Threshold = ", thresh), line = 1)
      
      ### subtitle with statistics
      title(sub =
              paste0("Prop. Sampled cells retained = ",
                     round(
                       land[land$Threshold == thresh, "SampledIncluded"] / 
                         max(land[, "SampledIncluded"], na.rm = TRUE),
                       3),
                     "\n",
                     "Prop. Occurrences excluded = ", 
                     round(
                       land[land$Threshold == thresh, "OccurrencesExcluded"],
                       2),
                     "\n",
                     "Sampled cells excuded = ", 
                     land[land$Threshold == thresh, "SampledExcluded"],
                     "\n", 
                     "Unsampled cells added = ", 
                     land[land$Threshold == thresh, "UnsampledAdded"]),
            line = 3.5)
      
      ### thresholded map
      terra::plot(atlas_boundary_thresh, 
                  colNA = rgb(0.5, 0.5, 0.5),
                  col = c(rgb(0.5, 0.5, 0.5), rgb(1, 0, 0)),
                  alpha = 0.35,
                  add = TRUE,
                  legend = FALSE)
      
      ### add in boundary polygon of thresholded map
      if(!is.null(boundary_poly)) {
        if(length(boundary_poly) >  0) {
          terra::plot(boundary_poly, add = TRUE)
        }
        if(length(boundary_poly) == 0) {
          rect(xleft   = terra::ext(atlas_raster_extend)$xmin,
               ybottom = terra::ext(atlas_raster_extend)$ymin,
               xright  = terra::ext(atlas_raster_extend)$xmax,
               ytop    = terra::ext(atlas_raster_extend)$ymax)
        }
      }
    }
  }
  
  output <- list(Thresholds = data.frame(All_Sampled      = 0,
                                         All_Occurrences  = occ_thresh,
                                         Gain_Equals_Loss = gain_loss_thresh,
                                         Sampled_Only     = 1),
                 Data =  land)
  return(output)
}
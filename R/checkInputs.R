################################################################################
# 
# inputChecking.R
# Version 2.0
# 18/05/2023
#
# Updates:
#   18/05/2023: function added (moved from upgrain and upgrain.threshold)
#               uses inherits for class conditions
#
# Error checking for inputs of upgrain and upgrain.threshold
#
################################################################################

checkInputs <- function(inputFunction,
                        atlas.data,
                        cell.width,
                        scales,
                        threshold,
                        thresholds,
                        method,
                        occupancies,
                        model,
                        models,
                        extent) {
  
  ##############################################################################
  ### both upgrain and upgrain.threshold
  
  if(inputFunction %in% c("upgrain.threshold", "upgrain")) {
    ### Error checking: scales
    if(scales < 2){
      stop("Scales must be >1 (at least three grain sizes needed for downscaling)",
           call. = FALSE)
    }
    
    ### Error checking: if data frame or sfc needs cell width
    if(inherits(atlas.data, "data.frame") | inherits(atlas.data, "sf")) {
      if(is.null(cell.width)) {
        stop("If data is data.frame cell.width is required", call. = FALSE)
      }
    }
    
    ### Error checking: if data frame check column names are correct
    if(inherits(atlas.data, "data.frame") & !inherits(atlas.data, "sf")) {
      if(ncol(atlas.data) != 3) {
        stop("Input data frame must contain three columns named 'x', 
           'y', and 'presence' in that order",
             call. = FALSE)
      }
      
      if(sum(names(atlas.data) != c("x", "y", "presence")) > 0) {
        stop("Input data frame must contain three columns named 'x', 
            'y', and 'presence' in that order",
             call. = FALSE)
      }
    }
    
    # ### Error checking: raster and sp packages no longer supported
    # if(inherits(atlas.data, "SpatialPointsDataFrame")) {
    #   stop('SpatialPointsDataFrame and sp package no longer supported as of 5.0-0:
    #      Convert to an sf object using st_as_sf(your_object))',
    #        call. = FALSE)
    # }
    # 
    # if(inherits(atlas.data, "RasterLayer")) {
    #   stop("'raster' package no longer supported as of 5.0-0: 
    #      Convert to a SpatRaster object using the 'terra' package instead",
    #        call. = FALSE)
    # }
  }
  
  ##############################################################################
  ### upgrain.threshold specific
  
  if(inputFunction == "upgrain.threshold") {
    ### Error checking: thresholds are between and include 0 and 1
    if(min(thresholds) != 0){
      stop("Threshold values must be between and include 0 and 1",
           call. = FALSE)
    }
    if(max(thresholds) != 1){
      stop("Threshold values must be between and include 0 and 1",
           call. = FALSE)
    }
  }
  
  ##############################################################################
  ### upgrain specific
  
  if(inputFunction == "upgrain") {
    ### Error checking: either a threshold is given or a method selected, not both
    if((!is.null(threshold)) & (!is.null(method))) {
      stop("Both a threshold and a model have been given. Please set one to NULL",
           call. = FALSE)
    }
    
    if(is.null(method)) {
      ### Error checking: threshold is between 0 and 1
      if(min(threshold) < 0){
        stop("Threshold value must be between 0 and 1", call. = FALSE)
      }
      if(max(threshold) > 1){
        stop("Threshold value must be between 0 and 1", call. = FALSE)
      }
    }
    
    if(is.null(threshold)) {
      ### Error checking: more than one method given
      if(length(method) > 1) {
        stop("More than one method selected", call. = FALSE)
      }
      
      ### Error checking: more than one method given
      if(sum(method %in% c("Sampled_Only",
                           "All_Sampled",
                           "All_Occurrences",
                           "Gain_Equals_Loss")) > 1) {
        stop("More than one method selected", call. = FALSE)
      }
      
      ### Error checking: method is one of the selections
      if(!method %in% c("Sampled_Only",
                        "All_Sampled",
                        "All_Occurrences",
                        "Gain_Equals_Loss")) {
        stop("Method is not one of the options", call. = FALSE)
      }
    }
  }
  
  ##############################################################################
  ### downscale
  
  if(inputFunction == "downscale") {
    if(inherits(occupancies, "upgrain")) {
      ### Error checking: input data frame correct
      if(ncol(occupancies) != 2) {
        stop("Input data must be a data frame with two columns (cell area and 
           occupancy", call. = FALSE)
      }
      
      ### Error checking: extent required
      if(is.null(extent)) {
        stop("Total extent required", call. = FALSE)
      }
      
      ### Error checking: extent larger than largest grain size
      if(extent < max(occupancies[, 2])) {
        stop("Total extent is smaller than the largest grain size! 
           Are the units correct?", call. = FALSE)
      }
      
      ### Error checking: occupancies are between 0 and 1
      if(min(occupancies[, 2]) < 0) {
        stop("Occupancies must be proportion of cells occupied (values must be
           between 0 - 1)", call. = FALSE)
      }
      if(max(occupancies[, 2]) > 1) {
        stop("Occupancies must be proportion of cells occupied (values must be
           between 0 - 1)", call. = FALSE)
      }
    }
    
    ### Error checking: model name is correct
    if (!model %in% c("Nachman", "PL", "Logis", "Poisson", "NB", "GNB", "INB",
                      "FNB", "Thomas")) {
      stop("Model name invalid", call. = FALSE)
    }
  }
  
  ##############################################################################
  ### hui.downscale
  
  if(inputFunction == "hui.downscale") {
    ### Error checking: if data frame requires extent and cell width
    if(inherits(atlas.data, "data.frame") & !inherits(atlas.data, "sf")) {
      if(is.null(extent)) {
        stop("Extent required if data input is data frame of coordinates",
             call. = FALSE)
      }
      if(is.null(cell.width)) {
        stop("If data is data.frame cell.width is required", call. = FALSE)
      }
    }
    
    ### Error checking: if sf requires extent and cell width
    if(inherits(atlas.data, "sf")) {
      if(is.null(extent)) {
        stop("Extent required if data input is sf spatial points object",
             call. = FALSE)
      }
      if(is.null(cell.width)) {
        stop("If data is sf spatial points object cell.width is required",
             call. = FALSE)
      }
    }
  }
  
  ##############################################################################
  ### ensemble.downscale
  
  if(inputFunction == "ensemble.downscale") {
    ### Error checking: at least one model selected
    if (length(models) == 1) {
      if(models != "all") {
        stop("Only one model selected: ensemble modelling not applicable", 
             call. = FALSE)
      }
    }
    
    ### Error checking: model names
    if(any(!models %in% c("Nachman", "PL", "Logis", "Poisson", "NB", "GNB", 
                          "INB", "FNB", "Thomas", "Hui", "all"))) {
      stop("One or more model names invalid", call. = FALSE)
    }
    
    ### If not an upgrain object
    if(!inherits(occupancies, "upgrain")) {
      ### Error checking: extent given
      if(is.null(extent)) {
        stop("No extent given and occupancies is not of class 'upgrain'",
             call. = FALSE)
      }
      
      ### Error checking: extent larger than largest grain size
      if(extent < max(occupancies[, 2])) {
        stop("Total extent is smaller than the largest grain size! 
           Are the units correct?", call. = FALSE)
      }
      
      ### Error checking: for hui model occupancies must be upgrain object
      if(length(models) == 1) {
        if(models == "all") {
          stop("Hui model can not be run if occupancies is not of class 'upgrain'",
               call. = FALSE)
        }
      }
      if(length(models) > 1) {
        if(any(models == "Hui")) {
          stop("Hui model can not be run if occupancies is not of class 'upgrain'",
               call. = FALSE)
        }
      }
    }
    
    ### Error checking: if model = Hui then cell.width and extent must be present
    if(inherits(occupancies, "upgrain")) {
      cell.width <- terra::res(occupancies$atlas.raster.stand)[1]
      extent <- occupancies$extent.stand
      
      if(length(models) == 1) {
        if(models == "all") {
          if(is.null(cell.width)) {
            stop("Cell.width must be specified for the Hui model", call. = FALSE)
          }
          if(is.null(extent)) {
            stop("Extent must be specified for the Hui model", call. = FALSE)
          }
        }
      }
      
      if(length(models) >  1) {
        if(any(models == "Hui")) {
          if(is.null(cell.width)) {
            stop("Cell.width must be specified for the Hui model", call. = FALSE)
          }
          if(is.null(extent)) {
            stop("Extent must be specified for the Hui model", call. = FALSE)
          }
        }
      }
    }
  }
}
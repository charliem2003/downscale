################################################################################
# 
# ensemble.downscale.R
# Version 2.0
# 22/05/2023
#
# Updates:
#   22/05/2023: v2.0 - CONVERTED TO TERRA AND SF
#               error checking moved to checkInputs
#               uses inherits for class conditions
#   16/05/2023: Simple reformatting
#   26/10/2021: uses on.exit to return to original par settings
#   24/11/2016: Bug fixed with apply in error checking
#   13/03/2015: if 0's predicted don't plot them
#   24/02/2015: means calculated as mean of log occupancies
#   24/02/2015: warning messages for inconsistent results
#   23/02/2015: different tolerances allowed for modelling and predicting
#   05/02/2015: plot = TRUE argument added
#   05/02/2015: improved warning control
#
# Ensemble modelling
#
# Args:
#   occupancies: data frame of observed area of occupancies and cell areas, or 
#                object from upgrain
#   new.areas: vector of grain sizes (in same units as occupancy) for which area 
#               of occupancy will be predicted.
#   extent: total area in same units as occupancy
#   tolerance_mod: The tolerance
#           used during integration in the Thomas model during optimisation of
#           parameters.
#   tolerance_pred: The tolerance used during the prediction stage. 
#   tolerance_hui: The tolerance used in the Hui model
#   models: vector of chosen downscaling models. Default models = "all" runs all
#           available models.
#   plot: if TRUE predictions of all models are plotted against grain size along
#         with the mean of all models.
#   verbose: if TRUE prints updates on modelling status.
# Returns:
#    a dataframe. The first column cell.area is the grain sizes used for 
#    predictions. The final column Means are the mean predictions of all models
#    for each grain size. Intermediate columns are the predicted occupancies for
#    the selected downscaling models.
#
################################################################################

ensemble.downscale <- function(occupancies,
                               new.areas,
                               extent = NULL,
                               cell.width = NULL,
                               models = "all",
                               tolerance_mod = 1e-6,
                               tolerance_pred = 1e-6,
                               tolerance_hui = 1e-6,
                               starting_params = NULL,
                               plot = TRUE,
                               verbose = TRUE) {

  ##############################################################################
  ### Error checking
  
  checkInputs(inputFunction = "ensemble.downscale",
              occupancies   = occupancies,
              models        = models,
              extent        = extent)

  ##############################################################################
  ### data handling
  
  ### list of models to run
  if(length(models) == 1) {
    if(models == "all") {
      model.list <- c("Nachman","PL","Logis","Poisson","NB",
                      "GNB","INB","FNB","Thomas", "Hui")
    }
  } else {
    model.list <- models    
  }
  
  ### extract info from upgrain object
  if(inherits(occupancies, "upgrain")) {
    cell.width <- terra::res(occupancies$atlas.raster.stand)[1]
    extent <- occupancies$extent.stand
  }
  
  ### any optional parameters for specific models
  starting_params_opts <- starting_params
  starting_params_mods <- NULL
  if(!is.null(starting_params)) {
    starting_params_mods <- names(starting_params_opts)
  }
  
  ##############################################################################
  ### Results storage
  
  all.predicted <- as.data.frame(matrix(NA, 
                                        ncol = (length(model.list) + 1),
                                        nrow = length(new.areas)))
  colnames(all.predicted) <- c("Cell.area", model.list)
  all.predicted$Cell.area <- new.areas
  
  ##############################################################################
  ### Modelling
  
  for (i in 1:length(model.list)) {
    model.run <- model.list[i]
    
    ### see if there are user-inputted starting parameters
    if(!is.null(starting_params_mods)) {
      if(any(starting_params_mods == model.run)) {
        starting_params_mod <- starting_params_opts[[model.run]]
      } else {
        starting_params_mod <- NULL
      }
    } else {
      starting_params_mod <- NULL
    }
    
    if(verbose == TRUE){
      cat(paste(model.run, "model is running..."))
    }
    
    ### For all models except Hui
    if(model.run != "Hui") {
      mod <- downscale(occupancies     = occupancies,
                       model           = model.run,
                       extent          = extent,
                       tolerance       = tolerance_mod,
                       starting_params = starting_params_mod)
      
      est <- predict.downscale(object    = mod,
                               new.areas = new.areas,
                               extent    = extent,
                               tolerance = tolerance_pred,
                               plot      = FALSE)
      all.predicted[, i + 1] <- est$predicted[, "Occupancy"]
    }
    
    ### For the Hui model
    if(model.run == "Hui") {
      new.areas.hui <- new.areas[new.areas < (cell.width ^ 2)]
      est <- hui.downscale(atlas.data = occupancies, 
                           cell.width = cell.width, 
                           new.areas  = new.areas.hui, 
                           extent     = extent,
                           tolerance  = tolerance_hui)
      all.predicted[1:length(new.areas.hui), i+1] <- est$predicted[,"Occupancy"]
    }
    
    if(verbose == TRUE){
      cat(paste("  complete", "\n"))
    }
  }
  
  ### save results
  all.predicted$Means <- exp(rowMeans(log(all.predicted[, -1]), na.rm = TRUE))
  aoo.predicted <- all.predicted
  aoo.predicted[, -1] <- aoo.predicted[, -1] * extent
  
  ##############################################################################
  ### optional plotting
  
  if (plot == TRUE) {
    parOrig <- par(no.readonly = TRUE)
    on.exit(par(parOrig))
    par(mfrow = c(3, ceiling(length(model.list) / 3)), mar = c(5, 5, 3, 1))
    
    for (i in 1:length(model.list)) {
      predicted <- all.predicted
      predicted[predicted == 0] <- NA
      plot(predicted[, 2] ~ all.predicted$Cell.area,
           type = "n",
           log  = "xy",
           xlim = c(min(c(all.predicted$Cell.area,
                          est$observed[, "Cell.area"]), na.rm = TRUE),
                    max(c(all.predicted$Cell.area,
                          est$observed[, "Cell.area"]), na.rm = TRUE)),
           ylim = c(min(predicted[, -1], na.rm = TRUE), 1),
           xlab = "Log cell area",
           ylab = "Log occupancy",
           main = paste(model.list[i], "model"))
      
      points(est$observed[, "Occupancy"] ~ est$observed[, "Cell.area"],
             type="b",
             lwd=2)
      points(predicted$Means ~ all.predicted$Cell.area,
             type="b",
             lwd=2,
             col = "dark grey")
      points(predicted[, i + 1] ~ all.predicted$Cell.area,
             type="b",
             lwd=2,
             col = "red")
    }
  }
  
  ##############################################################################
  ### Output
  
  output <- list(Occupancy = all.predicted,
                 AOO = aoo.predicted)
  return(output)
}


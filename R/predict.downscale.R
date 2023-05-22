################################################################################
# 
# predict.downscale.R
# Version 2.0
# 22/05/2023
#
# Updates:
#   16/05/2023: Simple reformatting
#               inherits for checking class conditions
#   26/10/2021: uses on.exit to return to original par settings
#   24/11/2016: Bug fixed with apply in error checking
#   13/03/2015: if 0's predicted don't plot them
#   03/02/2015: plot function added
#   03/02/2015: output defined as class 'downscale'
#   03/02/2015: observed data included in output
#   02/02/2015: Help file added to
#   30/01/2015: Thomas model added
#
# Predict area of occupancy at fine grain sizes using parameters in a downscale
# object estimated from coarse grain sizes using downscale
#
# Args:
#   mod.fit: model output of class 'downscale' (created from function downscale)
#   new.areas: vector of grain sizes for model prediction
#   extent: total area (same units as newdata)- required only for FNB and Thomas
#           models
#   tolerance: tolerance for integration of Thomas model. Lower numbers allow
#              for greater accuracy but require longer processing times
#   plot: if TRUE plots observed and predicted occupancies against grain size on
#         a log-log plot
# 
# Returns:
#   list of three objects of class 'predict.downscale'
#     $model      Downscaling model used
#     $predicted  Data frame of grain sizes and predicted occupancies
#     $observed   Data frame of grain sizes and observed occupancies used for 
#                 modelling
#
################################################################################

predict.downscale <- function(object, 
                              new.areas, 
                              tolerance = 1e-6, 
                              plot = TRUE,
                              ...) {
  mod.fit <- object
  
  ### error checking of inputs
  if(!inherits(mod.fit, "downscale")) {
    stop("Input data not of class 'downscale'")
  }
  
  ### extract parameters from 'downscale' object
  params <- as.list(mod.fit$pars)
  model <- mod.fit$model
  extent <- mod.fit$extent
  predict.function <- getFunction(paste("Predict", model, sep = ""))
  
  ### run predict functions for each model to estimate AOO
  if (model %in% c("Nachman", "PL", "Logis", "Poisson", "NB", "GNB", "INB")) {    
    AOO <- exp(predict.function(par  = params,
                                area = new.areas))
  }
  
  if (model == "FNB") {
    AOO <- exp(predict.function(par    = params, 
                                area   = new.areas, 
                                extent = extent))
  }
  
  if (model == "Thomas") {
    AOO <- exp(predict.function(par       = params,
                                tolerance = tolerance,
                                area      = new.areas, 
                                extent    = extent))
  }
  
  ### results data frame
  expected <- data.frame(Cell.area = new.areas,
                         Occupancy = AOO,
                         AOO       = AOO * extent)
  
  ### error checking in results
  if(anyNA(expected$Occupancy)) {
    warning("Predicted results may be innaccurate: one or more NA's predicted.")
  }
  
  if(sum(expected$Occupancy == 0, na.rm = TRUE) > 0) {
    warning("Predicted results may be innaccurate: one or more 0's predicted.")
  }
  
  if(!anyNA(expected$Occupancy) & nrow(expected) > 1) {
    if(any(sapply(expected$AOO, FUN = function(x) { AOO[x] > AOO[x + 1] }),
           na.rm = TRUE)) {
      warning("Scaling is inconsistent:
               \nlarger occupancies predicted at finer grain sizes.
               \nIf model = Thomas try a smaller tolerance value (e.g. 1e-7)")
    }
  }
  
  ### output results as class 'predict.downscale'
  output <- list(model     = model,
                 predicted = expected,
                 observed  = mod.fit$observed)
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

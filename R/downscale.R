################################################################################
# 
# downscale.R
# Version 2.0
# 19/05/2023
#
# Updates:
#   19/05/2023: v2.0 - CONVERTED TO TERRA AND SF
#               error checking moved to checkInputs
#               uses inherits for class conditions
#   16/05/2023: renaming of optimisePars functions. Simple reformatting
#   22/11/2021: basic reformatting
#   08/05/2015: extent now required
#   03/02/2015: output defined as class 'downscale'
#   03/02/2015: observed data included in output
#   02/02/2015: help file updated
#   02/02/2015: error check for model name and extent
#   30/01/2015: Thomas model added 
#
# Model area of occupancy against grain size for downscaling. 
#
# Args:
#   occupancies: data frame of observed area of occupancies and cell areas, or 
#                object from upgrain
#   model: function to use
#   extent: total area (same units as area)
#   tolerance: tolerance for integration of Thomas model. Lower numbers allow 
#              for greater accuracy but require longer processing times
#   starting_params: optional list of parameter values
#
# Returns:
#   list of three objects of class 'downscale'
#     $model    Downscaling model used
#     $pars     List of parameters estimated from optimisation procedure
#     $observed Dataframe of grain sizes and observed occupancies used for 
#               modelling
#
################################################################################

downscale <- function(occupancies,
                      model,
                      extent = NULL,
                      tolerance = 1e-6,
                      starting_params = NULL) {
  
  if(inherits(occupancies, "upgrain")) {
    extent <- occupancies$extent.stand
    occupancies <- occupancies$occupancy.stand[, -2]
  }
  
  ##############################################################################
  ### Error checking
  checkInputs(inputFunction = "downscale",
              occupancies = occupancies,
              model = model,
              extent = extent)
  
  ##############################################################################
  ### data manipulation
  input.data <- DataInput(occupancy = occupancies[, 2],
                          area      = occupancies[, 1],
                          extent    = extent)
  model <- model
  if(is.null(starting_params)) {
    starting_params <- NULL
  }
  
  ##############################################################################
  ### optimisation
  
  if(model %in% c("Nachman", "PL", "Logis", "Poisson", "NB", "INB")) {
    optim.pars <- suppressWarnings(
      OptimisePars(area = input.data$Cell.area[!is.na(input.data$Occ)],
                   observed = input.data$Occ[!is.na(input.data$Occ)],
                   model = model,
                   starting.params = starting_params))
  }
  
  if(model == "Logis") { 
    optim.pars <- suppressWarnings(
      OptimiseParsLogis(area = input.data$Cell.area[!is.na(input.data$Occ)], 
                        observed = input.data$Occ[!is.na(input.data$Occ)],
                        model = model,
                        starting.params = starting_params))
  }
  
  if(model == "GNB") { 
    optim.pars <- suppressWarnings(
      OptimiseParsGNB(area = input.data$Cell.area[!is.na(input.data$Occ)], 
                      observed = input.data$Occ[!is.na(input.data$Occ)],
                      model = model,
                      starting.params = starting_params))
  }
  
  if (model == "FNB") { 
    optim.pars <- suppressWarnings(
      OptimiseParsFNB(area = input.data$Cell.area[!is.na(input.data$Occ)], 
                      observed = input.data$Occ[!is.na(input.data$Occ)],
                      extent = extent,
                      model = model,
                      starting.params = starting_params))
  }
  
  if (model == "Thomas") { 
    optim.pars <- suppressWarnings(
      OptimiseParsThomas(area = input.data$Cell.area[!is.na(input.data$Occ)], 
                         observed = input.data$Occ[!is.na(input.data$Occ)],
                         extent = extent,
                         model = model,
                         tolerance = tolerance,
                         starting.params = starting_params))
  }
  
  ### output
  observed <- data.frame(Cell.area = input.data$Cell.area,
                         Occupancy = input.data$Occ)
  output <- list(model    = model,
                 pars     = unlist(optim.pars),
                 observed = observed,
                 extent   = extent)
  class(output) <- "downscale"
  return(output)
}

################################################################################
# 
# OptimiseParametersFNB.R
# Version 1.1
# 18/08/2016
#
# Updates:
#   18/08/2016: names of starting parameters harmonised with paper
#
# Optimisation procedure of finding parameters that best fit observed data for
# the Finite Negative Binomial model
#
# Args:
#   Area: Vector of grain sizes for obsvered area of occupancies
#   Observed: Vector of observed area of occupancies
#   extent: Extent of study area (km2)
#   model = "FNB"
#
# Returns:
#   optim.pars: list of parameters estimated from optimisation procedure
#
################################################################################

OptimiseParametersFNB <- function(area,
                                  observed,
                                  extent,
                                  model = "FNB",
                                  starting.params = NULL) {
  # Retrive residual function, downscaling function and starting parameters
  # for model of choice
  resid.fun <- getFunction(paste("Resid", model, sep = ""))
  pred.fun <- getFunction(paste("Predict", model, sep = ""))  
  if(is.null(starting.params)) {
    starting.pars <- get(paste("Params", model, sep = ""))
  }
  if(!is.null(starting.params)) {
    starting.pars <- starting.params
  }
  
  # Optimisation procedure
  optimisation <- minpack.lm::nls.lm(par = starting.pars,
                                     fn = resid.fun,
                                     area = area,
                                     extent = extent,
                                     observed = log(observed),
                                     lower = c("N" = -Inf, "k" = 0),
                                     upper = c("N" = Inf, 
                                               "k" = min(area, na.rm = TRUE) * 
                                                 100),
                                     control = minpack.lm::nls.lm.control(
                                       maxiter = 1000))
  optim.pars <- as.list(coef(optimisation))
  return(optim.pars)
}

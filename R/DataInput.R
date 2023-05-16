################################################################################
# 
# DataInput.R
# Version 1.2
# 16/05/2023
#
# Updates:
#   16/05/2023: Simple reformatting
#   06/10/2015: Bug in scale of endemism fixed
#
# Formatting of input data
#
################################################################################

DataInput <- function(occupancy,
                      area,
                      extent) {
  
  ### error checking: Area same length as occupancy
  if (length(occupancy) != length(area)) {
    stop("Number of grain sizes does not equal number of observed occupancies",
         call. = FALSE)
  }
  
  observed <- data.frame (Cell.area = area,
                          Occ = occupancy)
  observed <- observed[order(observed$Cell.area), ]
  
  ### Convert cell areas to proportion of total area
  ### (should standardise any data to fit within starting parameter range)
  #observed[, "Cell.area"] <- observed[, "Cell.area"] ^ 2 / extent
  
  ### Scale of saturation
  saturation.point <- suppressWarnings(min(which(observed$Occ == 1)))
  if (saturation.point < length(observed$Occ)) {
    observed[(saturation.point + 1):length(observed$Occ), "Occ"] <- NA
  }
  
  ### error checking: Enough scales before scale of saturation
  if ((saturation.point == 1) | (saturation.point == 2)) {
    stop("Not enough scales before scale of saturation for modelling",
         call. = FALSE)
  }
  
  ### Scale of endemism
  no_cells <- extent / observed$Cell.area
  endemism.point <- suppressWarnings(
    min(which(observed$Occ == (1 / no_cells))))
  if (endemism.point < length(observed$Occ)) {
    observed[(endemism.point + 1):length(observed$Occ), "Occ"] <- NA
  }

  ### error checking: Enough scales before scale of endemism
  if ((endemism.point == 1) | (endemism.point == 2)) {
    stop("Not enough scales before scale of endemism for modelling",
         call. = FALSE)
  }
  return(observed)
}


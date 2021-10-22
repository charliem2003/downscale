\name{hui.downscale}
\alias{hui.downscale}
\title{Predict occupancy at fine grain sizes using the Hui model}
\usage{
hui.downscale(atlas.data, cell.width, new.areas, extent = NULL,
              tolerance = 1e-6, plot = FALSE)
}
\arguments{
\item{atlas.data}{either an object of class \code{"upgrain"}; or a raster file where 1 = presence and 0 = absence; or an object of class SpatialPointsDataFrame atlas data where the data frame contains a single column of presence (1) and absences (0); or a data frame containing these columns: \code{presence} = presence - absence data; \code{x} = easting coordinates of cells; \code{y} = northing coordinates of cells.}

\item{cell.width}{the cell width of the atlas data.}

\item{new.areas}{vector of grain sizes as cell area (in same units as cell.width but squared e.g. km^{2}) at fine scales for model prediction.}

\item{extent}{the extent of the atlas data in same units as \code{cell.areas}. If the input data is of class \code{"upgrain"} this can be left as \code{NULL}.}

\item{tolerance}{tolerance for root solving to estimate probability of absence at the fine scale (default = \code{1e-6}).}

\item{plot}{if \code{plot = TRUE} (default = \code{FALSE}) plots observed and predicted occupancies against grain size on a log-log plot.}
}
\value{
Returns an object of class \code{'predict.downscale'} with three
  objects:
   \item{model}{Downscaling model used (Hui model in this case).}
   \item{predicted}{Data frame containing two columns:
   \tabular{lll}{
     \code{Cell.area} \tab  \tab Grain sizes for which occupancy have been estimated\cr
     \code{Occupancy} \tab  \tab Predicted area of occupancy for each grain size\cr
     }
   }
   \item{observed}{Data frame containing two columns:
   \tabular{lll}{
     \code{Cell.area} \tab  \tab Grain sizes for which occupancy have been
       observed\cr
     \code{Occupancy} \tab  \tab Observed area of occupancy for each grain
       size\cr
     }
   }
}
\description{
Predict area of occupancy at fine grain sizes from atlas data
  using the Hui model. There is also a simple plot function. Essentially the function is equivalent to \code{\link{downscale}} and \code{\link{predict.downscale}} combined.
}
\details{
The function estimates fine-scale occupancy from atlas-scale data
  using the model developed by Hui. The model requires presence-absence data
  only at a single scale (the atlas scale) and calculates the probability of
  occupancy at a fine scale through the conditional probability that at the
  atlas scale a randomly chosen cell adjacent to an occupied cell is also
  occupied. Plotting can be called directly from
  \code{\link{hui.downscale}} or from \code{\link{plot.predict.downscale}}.
  
NOTE: if comparing occupancy predictions from the Hui model with those from other models using \code{\link{downscale}} it is important that the atlas data used here is the standardised atlas data whose extent is the same as the largest grain size used in the other models. This ensures that all models are modelling the proportion of occupancy of the same extent. This standardised atlas raster can be obtained through the \code{\link{upgrain}} function (see \code{vignette("Upgraining", package = "downscale")} for more details.
}
\examples{
## example species data
data.file <- system.file("extdata", "atlas_data.txt", package = "downscale")
atlas.data <- read.table(data.file, header = TRUE)

## if the input data is a data frame it must have the columns "x", "y"
## and "presence"
head(atlas.data)

## Fit Hui model to atlas.data
(hui <- hui.downscale(atlas.data, 
                      cell.width = 10,
                      extent = 228900,
                      new.areas = c(1, 2, 5, 15, 50,75),
                      plot = TRUE))

## Fit Hui model to standardised atlas data for comparison with other models
## First, upgrain data (using All Sampled threshold)
occupancy <- upgrain(atlas.data,
                     cell.width = 10,
                     scales = 3,
                     method = "All_Sampled")
                     
## the "upgrain" object can be used as input for the Hui model
(hui.stand <- hui.downscale(occupancy,
                            cell.width = 10,
                            new.areas = c(1, 2, 5, 15, 50,75),
                            plot = TRUE))

## compare the area of occupancy (AOO) predictions of the two models
hui$predicted
hui.stand$predicted
}
\author{
Charles Marsh <\email{charliem2003@gmail.com}> with input from
  Louise Barwell and Cang Hui.
}
\references{
Hui, C. (2009). On the scaling patterns of species spatial
  distribution and association. \emph{Journal of Theoretical Biology} 261,
  481-487.

Hui, C., McGeoch, M.A., & Warren, M. (2006). A spatially explicit
  approach to estimating species occupancy and spatial correlation.
  \emph{Journal of Animal Ecology} 7, 140-147.
  
  Groom, Q., Marsh, C.J., Gavish, Y. Kunin, W.E. (2018). How to predict fine resolution occupancy from coarse occupancy data, \emph{Methods in Ecology and Evolution}. In press.
  
  Marsh, C.J, Barwell, L.J., Gavish, Y., Kunin, W.E. (2018). downscale: An R package for downscaling species occupancy from coarse-grain data to predict occupancy at fine-grain sizes, \emph{Journal of Statistical Software, Code Snippets} 86(3), 1-20.
}
\seealso{
See \code{\link{downscale}} for estimating parameters of a
  downscaling function from observed occupancies at coarse grain sizes using
  other downscaling models.
  
See \code{\link{upgrain}} for creating extent-standardised atlas data as an input.

See \code{\link{ensemble.downscale}} for ensemble modelling of multiple downscaling models.
}


\name{predict.downscale}
\alias{predict.downscale}
\title{Predict occupancy at fine grain sizes}
\usage{
\method{predict}{downscale}(object, new.areas, tolerance = 1e-6, plot = TRUE, ...)
}
\arguments{
\item{object}{a fitted object of class \code{'downscale'}.}

\item{new.areas}{vector of grain sizes (in squared units e.g. km^{2}) for which
area of occupancy will be predicted.}

\item{tolerance}{only applicable for the \code{Thomas} model. The tolerance
  used during integration in the Thomas model during optimisation of
  parameters. Lower numbers allow for greater accuracy but require longer
  processing times (default = \code{1e-6}).}

\item{plot}{if \code{plot = TRUE} (default) plots observed
  and predicted occupancies against grain size on a log-log plot using 
  \code{\link{plot.predict.downscale}}.}

\item{\dots}{arguments, including graphical parameters for \code{\link{plot.predict.downscale}}, passed to other methods.}
}
\value{
\code{predict} returns an object of class \code{'predict.downscale'} with three
  objects:
   \item{model}{Downscaling model used.}
   \item{predicted}{Data frame containing two columns:
   \tabular{lll}{
     \code{Cell.area} \tab  \tab Grain sizes for which occupancies have been
       estimated\cr
     \code{Occupancy} \tab  \tab Predicted proportion of occupancy for each grain
       size\cr
     \code{AOO} \tab  \tab Predicted area of occupancy (proportion of area occupancy multiplied by extent)\cr
     }
   }
   \item{observed}{Data frame containing two columns:
   \tabular{lll}{
     \code{Cell.area} \tab  \tab Grain sizes for which occupancies have been
       observed\cr
     \code{Occupancy} \tab  \tab Observed proportion of occupancy for each grain
       size\cr
     }
   }
}
\description{
Predict proportion of occupancy at fine grain sizes using parameters from an object of class \code{downscale} estimated from coarse grain sizes using \code{\link{downscale}}. Proportion of occupancy is converted to area of occupancy (AOO) by multiplying by the extent. There is also a simple plot function.
}
\details{
The function takes the parameters for a downscaling model estimated
  through \code{\link{downscale}} and uses the model to predict area of
  occupancy at finer grain sizes. See \code{\link{downscale}} for details on
  the downscaling models and their parameterisation. Plotting can be called
  directly from \code{\link{predict.downscale}} or from \code{\link{plot.predict.downscale}}.
  
  For predictions using the Thomas model, if the tolerance value is not sufficiently low it may lead to innacurate results. Typically, this will be indicated by fine grain sizes with higher predicted occupancies than those at larger grain sizes. In these cases try a lower tolerance value in the arguments.
}
\examples{
## example species data
data.file <- system.file("extdata", "atlas_data.txt", package = "downscale")
atlas.data <- read.table(data.file, header = TRUE)

## if the input data is a data frame it must have the columns "lon", "lat"
## and "presence"
head(atlas.data)

## explore thresholds using upgrain.threshold
thresh <- upgrain.threshold(atlas.data = atlas.data,
                            cell.width = 10,
                            scales = 3,
                            thresholds = seq(0, 1, 0.1))

## upgrain data (using All Sampled threshold)
occupancy <- upgrain(atlas.data,
                     cell.width = 10,
                     scales = 3,
                     method = "All_Sampled")

## Logistic model
(logis <- downscale(occupancies = occupancy,
                      model = "Logis"))

## Predictions occupancy at finer scales
logis.pred <- predict(logis,
                      new.areas = c(1, 2, 5, 25, 100),
                      plot = FALSE)
logis.pred
plot(logis.pred)

## It can also be a good idea to include the observed cell sizes in the
## predictions to observe the fit of the model (this time run plot = TRUE)
predict(logis,
        new.areas = c(1, 2, 5, 25, 100, 400, 1600, 6400),
        plot = TRUE)
}
\author{
Charles Marsh <\email{charliem2003@gmail.com}> with input from
  Louise Barwell.
}
\references{
Azaele, S., Cornell, S.J., & Kunin, W.E. (2012). Downscaling
  species occupancy from coarse spatial scales. \emph{Ecological
  Applications} 22, 1004-1014.

Barwell, L.J., Azaele, S., Kunin, W.E., & Isaac, N.J.B. (2014).
  Can coarse-grain patterns in insect atlas data predict local occupancy?
  \emph{Diversity and Distributions} 20, 895-907.
}
\seealso{
See \code{\link{downscale}} for estimating parameters of a
  downscaling function from observed occupancies at coarse grain sizes.
}


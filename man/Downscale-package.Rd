\name{downscale-package}
\alias{downscale-package}
\docType{package}
\title{
Downscaling Species Occupancy
}
\details{
\tabular{ll}{
Package: \tab downscale\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-05-18\cr
License: \tab GPL-2\cr
}
}
\description{
The package provides a set of functions that model the occupancy-area relationship (OAR) of known coarse scale data. The models are then extrapolated to predict the proportion of occupied area at finer grain sizes.
}
\section{Overview}{
The package provides three sets of functions for each stage of analysis:

1) \code{\link{upgrain}} and \code{\link{upgrain.threshold}} prepare atlas data for downscaling. 

2) \code{\link{downscale}} and \code{\link{hui.downscale}} model the OAR to the prepared data for one of ten possible downscaling models. 

3) \code{\link{predict.downscale}} and \code{\link{plot.predict.downscale}} take the model outputs and predict occupancy at finer grains. 

Finally, \code{\link{ensemble.downscale}} will run \code{\link{downscale}} and \code{\link{predict.downscale}} for a number of selected downscaling functions and calculate the mean predicted occupancies across all models.

The general flow of the package, and the inputs required for each function, is as follows:

\if{html}{\figure{Flow.png}{options: width="65\%" alt="Figure: Flow.png"}}
\if{latex}{\figure{Flow.png}{options: width=10cm}}

Two vignettes are available to guide users. Both work through examples in code:\cr

\code{vignette("Downscaling", package = "downscale")}\cr

\code{vignette("Upgraining", package = "downscale")}
}
\section{Credits}{
This package was created as part of deliverable D3.2 of WP3 of the project:
\bold{EU-BON: Building the European Biodiversity Observation Network} -
a 7th Framework Programme funded by the European Union under Contract No. 308454.
}

\author{
Charles Marsh with input from Louise Barwell and Cang Hui.

Maintainer: Charles Marsh <charliem2003@gmail.com>

For reporting bugs or requesting information please include \emph{'downscale'} in the subject line.
}
\references{
Azaele, S., Cornell, S.J., & Kunin, W.E. (2012). Downscaling species occupancy from coarse spatial scales. \emph{Ecological Applications} 22, 1004-1014.

Barwell, L.J., Azaele, S., Kunin, W.E., & Isaac, N.J.B. (2014). Can coarse-grain patterns in insect atlas data predict local occupancy? \emph{Diversity and Distributions} 20, 895-907.

Hui, C. (2009). On the scaling patterns of species spatial
  distribution and association. \emph{Journal of Theoretical Biology} 261,
  481-487.

Hui, C., McGeoch, M.A., & Warren, M. (2006). A spatially explicit
  approach to estimating species occupancy and spatial correlation.
  \emph{Journal of Animal Ecology} 7, 140-147.
  
Groom, Q., Marsh, C.J., Gavish, Y. Kunin, W.E. (2018). How to predict fine resolution occupancy from coarse occupancy data, \emph{Methods in Ecology and Evolution}. In press.

Marsh, C.J, Barwell, L.J., Gavish, Y., Kunin, W.E. (2018). downscale: An R package for downscaling species occupancy from coarse-grain data to predict occupancy at fine-grain sizes, \emph{Journal of Statistical Software, Code Snippets} 86(3), 1-20.
}
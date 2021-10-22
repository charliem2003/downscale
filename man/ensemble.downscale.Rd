\name{ensemble.downscale}
\alias{ensemble.downscale}
\title{Ensemble modelling of multiple downscaling functions}
\usage{
ensemble.downscale(occupancies, new.areas, extent, cell.width = NULL, 
                   models = "all", tolerance_mod = 1e-6, tolerance_pred = 1e-6,
                   tolerance_hui = 1e-6, starting_params = NULL, plot = TRUE, 
                   verbose = TRUE)
}
\arguments{
\item{occupancies}{Either a data frame containing two columns or an object of class "upgrain" from the upgrain function. If using the \code{Hui} model (or \code{model = "all"}) occupancies must be of class \code{"upgrain"}. If using a data frame the first column must be the grain sizes (cell area in squared units e.g. km^2). The second column is the proportion of occupancies at each grain size.}

\item{new.areas}{vector of grain sizes (in squared units e.g. km^{2}) for which
area of occupancy will be predicted.}

\item{extent}{total area in same units as occupancy. If using an object of class "upgrain", this is automatically inputted.}

\item{cell.width}{the cell width of the atlas data. Only required for \code{Hui} model.}

\item{tolerance_mod}{only applicable for the \code{Thomas} model. The
tolerance used during integration in the Thomas model during optimisation
of parameters. Lower numbers allow for greater accuracy but require longer
processing times (default = \code{1e-6}).}

\item{tolerance_pred}{only applicable for the \code{Thomas} model. The
tolerance used during the prediction stage.}

\item{tolerance_hui}{only applicable for the \code{Hui} model. The
tolerance used during integration in the Thomas model during optimisation
of parameters. Lower numbers allow for greater accuracy but require longer
processing times (default = \code{1e-6}).}

\item{models}{vector of chosen downscaling models. Default \code{models =
"all"} runs all available models. See \code{\link{downscale}} for list of
available models.}

\item{starting_params}{Starting values for model parameters if a model is not converging sufficiently. A list where each model to be specified is a list of parameter values. The parameter names must be the same for the default values (see \code{Details} for information on the parameters and how to specify them).}

\item{plot}{if \code{TRUE} predictions of all models (red) are plotted against
grain size along with the mean of all models (grey) and observed occupancies (black).}

\item{verbose}{if \code{TRUE} prints updates on modelling status.}
}
\value{
Returns a list of two data frames: 
\code{Occupancy} = proportion of occupancies; 
\code{AOO} = occupancies converted to area of occupancy. 

In each data frame the first column \code{cell.area} are the grain sizes used for predictions. The final column \code{Means} are the mean of the logged predictions of all models for each grain size. Intermediate columns are the predicted occupancies for the selected downscaling models.
}
\description{
Predict area of occupancy at fine grain sizes for multiple downscaling methods using \code{\link{downscale}} and \code{\link{predict.downscale}}. Occupancies are converted to area of occupancy (AOO) by multiplying by the total extent. The mean of the logged predicted occupancies of all models is then calculated.
}
\details{Ten downscaling models are available: \code{"Nachman"}, \code{"PL"}, \code{"Logis"}, \code{"Poisson"}, \code{"NB"}, \code{"GNB"}, \code{"INB"}, \code{"FNB"}, \code{"Thomas"} and \code{"Hui"}. They can be input in any order, or all models run through \code{models = "all"}. If the \code{Hui} model is included the input data must be an object of class \code{"upgrain"} generated through \code{\link{upgrain}}. See \code{\link{downscale}} and \code{\link{hui.downscale}} for more details of the available models.

The optimisation procedure requires initial starting values for all model parameters. In most cases the default values should work, however if one or more models are not converging adequately it is possible to input the starting parameters. For each model we wish to specify, the parameters must be in the form of a list with the same parameter names as in the table below. \code{starting_params} is then a list of these lists - the names of the lists are the same as the models. For example, if we wish to specify the starting parameters for the Nachman model and the Generalised Negative Binomial model the code would be:

\if{latex}{
\code{starting_params = list(Nachman = list("C" = 0.1, "z" = 0.01),} \cr
                       \code{GNB = list("C" = 0.1, "z" = 1, "k" = 0.01))}
}

\if{html}{
\verb{starting_params = list(Nachman = list("C" = 0.1, "z" = 0.01),
                       GNB = list("C" = 0.1, "z" = 1, "k" = 0.01))}
}

Please take particular note of captials. The default starting parameters are:

\if{html}{
 \tabular{lllllllllll}{
\code{"Nachman"} \tab \tab \tab C = 0.01; \tab \tab \tab z = 0.01 \tab \tab \tab \cr
\code{"PL"} \tab \tab \tab C = 0.01; \tab \tab \tab z = 0.01 \tab \tab \tab \cr
\code{"Logis"} \tab \tab \tab C = 0.01; \tab \tab \tab z = 0.01 \tab \tab \tab \cr
\code{"Poisson"} \tab \tab \tab gamma = 1e-8 \tab \tab \tab \tab \tab \tab \tab \cr
\code{"NB"} \tab \tab \tab C = 0.01; \tab \tab \tab k = 0.01 \tab \tab \tab \cr
\code{"GNB"} \tab \tab \tab C = 0.00001; \tab \tab \tab z = 1; \tab \tab k = 0.01\cr
\code{"INB"} \tab \tab \tab C = 1; \tab \tab \tab gamma = 0.01; \tab \tab b = 0.1\cr
\code{"FNB"} \tab \tab \tab N = 10; \tab \tab \tab k = 10 \tab \tab \tab \cr
\code{"Thomas"} \tab \tab \tab rho = 1e-8; \tab \tab \tab mu = 10; \tab \tab sigma = 1\cr
}
}
\if{latex}{\figure{Params.png}{options: width=7cm}}
}

\examples{
## example species data
data.file <- system.file("extdata", "atlas_data.txt", package = "downscale")
atlas.data <- read.table(data.file, header = TRUE)

## if the input data is a data frame it must have the columns "x", "y"
## and "presence"
head(atlas.data)

## upgrain data (using All Samples threshold)
occupancy <- upgrain(atlas.data,
                     cell.width = 10,
                     scales = 3,
                     method = "All_Sampled",
                     plot = FALSE)

## ensemble downscaling with an object of class upgrain
ensemble.downscale(occupancies = occupancy,
                   new.areas = c(1, 2, 5, 15, 50, 100, 400, 1600, 6400),
                   cell.width = 10,
                   models = c("Nachman", "PL", "Logis", "GNB", "FNB", "Hui"),
                   plot = TRUE)

## ensemble modelling with data frame of occupancies (not applicable for Hui
## model) with hypothetical species occupancy data
occupancy.dd <- data.frame(Cell.areas = c(100, 400, 1600, 6400),
                           Occupancy = c(0.16, 0.36, 0.59, 0.86))

## now extent must be specified (but cell.width not needed)
ensemble.downscale(occupancies = occupancy.dd,
                   new.areas = c(1, 2, 5, 15, 50, 100, 400, 1600, 6400),
                   extent = 384000,
                   models = c("Nachman", "PL", "Logis", "GNB", "FNB"),
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
  
Groom, Q., Marsh, C.J., Gavish, Y. Kunin, W.E. (2018). How to predict fine resolution occupancy from coarse occupancy data, \emph{Methods in Ecology and Evolution}. In press.
  
  Marsh, C.J, Barwell, L.J., Gavish, Y., Kunin, W.E. (2018). downscale: An R package for downscaling species occupancy from coarse-grain data to predict occupancy at fine-grain sizes, \emph{Journal of Statistical Software, Code Snippets} 86(3), 1-20.
}
\seealso{
See \code{\link{upgrain}} for creating objects of class \code{"upgrain"}. See \code{\link{downscale}} and \code{\link{predict.downscale}}, and \code{link{hui.downscale}} for downscaling models individually.
}


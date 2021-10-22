\name{downscale}
\alias{downscale}
\title{Model area of occupancy against grain size for downscaling}
\usage{
downscale(occupancies, model, extent = NULL, tolerance = 1e-06,
          starting_params = NULL)
}
\arguments{
\item{occupancies}{Either a data frame containing two columns or an object of class \code{"upgrain"} from the \code{\link{upgrain}} function. If using a data frame the first column must be the grain sizes (cell area in squared units e.g. km^{2}). The second column is the proportion of occupancies at each grain size.}

\item{model}{selected downscaling model, chosen from one of \code{"Nachman"},
\code{"PL"}, \code{"Logis"}, \code{"Poisson"}, \code{"NB"}, \code{"GNB"},
\code{"INB"}, \code{"FNB"}, \code{"Thomas"}. See \code{Details} below for
model descriptions.}

\item{extent}{total area in same units as occupancy. If using an object of class \code{"upgrain"}, this is automatically inputted.}

\item{tolerance}{only applicable for the \code{Thomas} model. The tolerance
  used during integration in the Thomas model during optimisation of
  parameters. Lower numbers allow for greater accuracy but require longer
  processing times (default = \code{1e-6}).}
  
\item{starting_params}{a list of starting values for model parameters. Useful if the default values are not converging during optimisation. The parameter names must be the same for the default values (see \code{Details} for information on the parameters).}
}

\value{
\code{downscale} returns an object of class \code{"downscale"} containing
  four objects:
  \item{model}{Downscaling model selected.}
  \item{pars}{Estimated parameters for the downscaling model.}
   \item{observed}{Data frame containing two columns:
   \tabular{lll}{
     \code{Cell.area} \tab  \tab Grain sizes for which occupancy have been
       observed.\cr
     \code{Occupancy} \tab  \tab Observed area of occupancy for each grain
       size.\cr
     }
   }
  \item{extent}{Only for \code{FNB} and \code{Thomas} models.}
}
\description{
Fits the log observed proportion of occupancies against log grain size (cell area) for coarse-scale data (typically atlas data) for nine commonly used downscaling models (see Azaele et al. 2012 and Barwell et al. 2014). See \code{\link{hui.downscale}} for downscaling using the Hui model. The parameters of the fitted models may then be used to estimate the area of occupancy at finer grain sizes than the observed data using \code{\link{predict.downscale}}. Presence-absence atlas data can be prepared for downscaling using \code{\link{upgrain}}.
}
\details{
Nine downscaling models are available. \code{area} is the grain size (cell area) and \code{extent} the total area in the same units:

\if{html}{
 \tabular{llll}{
   \code{"Nachman"} \tab  \tab Nachman model \tab
     \emph{log(1 - exp(-C * area ^ z))}\cr
   \code{"PL"} \tab  \tab  Power law model \tab
     \emph{log(C * area ^ z)}\cr
   \code{"Logis"} \tab  \tab  Logistic model \tab
     \emph{log((C * (area ^ z)) / (1 + (C * (area ^ z))))}\cr
   \code{"Poisson"} \tab  \tab  Poisson model \tab
     \emph{log(1 - (exp(-gamma * area)))}\cr
   \code{"NB"} \tab  \tab  Negative binomial model \tab
     \emph{log(1 - (1 + (gamma * area) / k) ^ -k)}\cr
   \code{"GNB"} \tab  \tab  Generalised negative binomial model \tab
     \emph{log(1 - (1 + (C * area ^ z) / k) ^ -k)}\cr
   \code{"INB"} \tab  \tab  Improved negative binomial model \tab
     \emph{log(1 - ((C * area ^ (b - 1)) ^ ((gamma * area) /
       (1 - C * area ^ (b - 1)))))}\cr
   \code{"FNB"} \tab  \tab  Finite negative binomial model \tab
     \emph{log(1 -
       ((gamma(N + ((extent * k) / area) - k) * 
       gamma(extent * k) / area) / 
       (gamma(N + ((extent * k) / area)) * 
       gamma(((extent * k) / area) - k)))
       }\cr
   \code{"Thomas"} \tab  \tab  Thomas model \tab see below\cr
   }
}
\if{latex}{\figure{Equations.png}{options: width=13cm}}

The finite negative binomial model (\code{"FNB"}) incorporates several gamma
functions. This may result in integers larger than is possible to store in
R. Therefore multiple precision floating point numbers (\code{\link{mpfr}}
function in package \pkg{Rmpfr}) are used to make calculations possible.

The Thomas model incorporates spatial point processes in order to estimate 
species aggregations. This involves multi-dimensional integration which may be 
time-consuming. Users can alter the tolerance value during the integration
process - a smaller value will give more accurate estimates but longer
processing times.
 
The optimisation procedure requires initial starting values for all model parameters. In most cases the default values should work, however if the model is not converging adequately it is possible to input the starting parameters. The parameters must be in the form of a list with the same parameter names as in the table below. For example for the \code{Nachman} model the code would be \code{starting_params = list("C" = 0.1, "z" = 0.01)}. Please take particular note of captials. The default starting parameters are:

\if{html}{
 \tabular{lllllllllll}{
\code{"Nachman"} \tab \tab \tab C = 0.01; \tab \tab \tab z = 0.01 \tab \tab \tab \cr
\code{"PL"} \tab \tab \tab C = 0.01; \tab \tab \tab z = 0.01 \tab \tab \tab \cr
\code{"Logis"} \tab \tab \tab C = 0.01; \tab \tab \tab z = 0.01 \tab \tab \tab \cr
\code{"Poisson"} \tab \tab \tab gamma = 1e-8 \tab \tab \tab \tab \tab \tab \tab \cr
\code{"NB"} \tab \tab \tab gamma = 0.01; \tab \tab \tab k = 0.01 \tab \tab \tab \cr
\code{"GNB"} \tab \tab \tab C = 0.00001; \tab \tab \tab z = 1; \tab \tab k = 0.01\cr
\code{"INB"} \tab \tab \tab C = 1; \tab \tab \tab gamma = 0.01; \tab \tab b = 0.1\cr
\code{"FNB"} \tab \tab \tab N = 10; \tab \tab \tab k = 10 \tab \tab \tab \cr
\code{"Thomas"} \tab \tab \tab rho = 1e-8; \tab \tab \tab mu = 10; \tab \tab sigma = 1\cr
}
}
\if{latex}{\figure{Params.png}{options: width=7cm}}

NOTE: for downscaling it is important that occupancies above the scale of saturation (the grain size at which all cells are occupied) and the scale of saturation (the grain size where only a single cell is occupied) are not included for modelling. The \code{downscale} functions will automatically set these occupancies to NA.
}
 
\examples{
## example species data
data.file <- system.file("extdata", "atlas_data.txt", package = "downscale")
atlas.data <- read.table(data.file, header = TRUE)

## if the input data is a data frame it must have the columns "x", "y"
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
                    
### predict occupancy at finer grain sizes
pred <- predict(logis,
                new.areas = c(1, 2, 5, 25, 100, 400, 1600, 6400))
                
### Plot predictions
plot(pred)

## Improved Negative Binomial model
(inb <- downscale(occupancies = occupancy,
                  model = "INB"))

## Specifying the starting parameters (gives a poorer fit in this case)
new.params <- list("C" = 0.1, "gamma" = 0.00001, "b" = 0.1)
(inb.new <- downscale(occupancies = occupancy,
                      model = "INB",
                      starting_params = new.params))

## plot the predictions of two FNB models using predict.downscale
predict(inb,
        new.areas = c(1, 2, 5, 25, 100, 400, 1600, 6400),
        plot = TRUE)
predict(inb.new,
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
  
Groom, Q., Marsh, C.J., Gavish, Y. Kunin, W.E. (2018). How to predict fine resolution occupancy from coarse occupancy data, \emph{Methods in Ecology and Evolution}. In press.

Marsh, C.J, Barwell, L.J., Gavish, Y., Kunin, W.E. (2018). downscale: An R package for downscaling species occupancy from coarse-grain data to predict occupancy at fine-grain sizes, \emph{Journal of Statistical Software, Code Snippets} 86(3), 1-20.
}
\seealso{
See \code{\link{upgrain}} for the preparation of presence-absence atlas data to occupancy data at several spatial scales.

The function output may be used as the input for \code{\link{predict.downscale}} for extrapolating downscaling functions to smaller grain sizes using the estimated parameters from the downscale output.

See \code{\link{hui.downscale}} for downscaling using the Hui model.
}

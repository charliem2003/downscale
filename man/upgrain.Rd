\name{upgrain}
\alias{upgrain}
\title{Upgraining of atlas data to larger grain sizes}
\usage{
upgrain(atlas.data, cell.width = NULL, scales, threshold = NULL,
        method = "All_Sampled", plot = TRUE, return.rasters = FALSE)
}
\arguments{
\item{atlas.data}{either a raster file of presence-absence; an object of class SpatialPointsDataFrame atlas data where the data frame contains a single column of presence (1) and absences (0); or a data frame of sampled cells, which must contain columns named 'x', 'y' and 'presence'.}

\item{cell.width}{if data is a data frame or SpatialPointsDataFrame, the cell widths of sampled cells. If data is a raster then leave as default (\code{= NULL})}

\item{scales}{the number of cells to upgrain. Upgraining will happen by
  factors of 2 - ie if scales = 3, the atlas data will be aggregated in 2x2
  cells, 4x4 cells and 8x8 cells.}

\item{threshold}{default = \code{NULL}. A user defined threshold for the proportion of unsampled atlas cells allowed within a cell at the largest grain size. Note: if a \code{method} is selected then threshold must be \code{NULL}.}

\item{method}{one of \code{"All_Sampled"}, \code{"All_Occurrences"}, \code{"Gain_Equals_Loss"} or \code{"Sampled_Only"} (default = \code{"All_Sampled"}). If the user wishes to define their own \code{threshold} then method must equal \code{NULL}. See \code{details} and \code{\link{upgrain.threshold}} for descriptions of the different methods.}

\item{plot}{default = \code{TRUE}. Plots the original atlas data alongside the standardised atlas data at each grain size.}

\item{return.rasters}{default = \code{FALSE}. If \code{TRUE} returns the extent-standardised atlas data upgrained to all grain sizes (NOTE: the extent-standardised atlas data at the original grain size is always returned regardless).}
}

\value{
Returns a list of class \code{"upgrain"} that can be used as a direct input to \code{\link{downscale}} and \code{\link{hui.downscale}}. The list contains five objects:
\item{threshold}{The proportion of unsampled atlas cells allowed within a cell at the largest grain size, either defined by the user through \code{threshold} or calculated according to \code{method}.}
\item{extent.stand}{The standardised extent after upgraining (equal to the extent of the largest grain size)}
\item{occupancy.stand}{Occupancy for each grain size where extent has been standardised. The data frame contains three columns:
   \tabular{lll}{
     \code{Cell.area} \tab  \tab Grain sizes for for each upgrained scale.\cr
     \code{Extent} \tab  \tab Extent for each grain size.\cr
     \code{Occupancy} \tab  \tab Observed area of occupancy for each grain
       size.\cr
     }
}
\item{occupancy.orig}{Original occupancies for each grain size before extent has been standardised. The data frame contains three columns:
   \tabular{lll}{
     \code{Cell.area} \tab  \tab Grain sizes for for each upgrained scale.\cr
     \code{Extent} \tab  \tab Extent for each grain size.\cr
     \code{Occupancy} \tab  \tab Observed area of occupancy for each grain
       size.\cr
     }
}
\item{atlas.raster.stand}{A raster layer of the extent-standardised atlas data}
\item{scaled.rasters}{If \code{return.rasters = TRUE}. A list containing the extent-standardised atlas data upgrained to all grain sizes}
}
\description{
Takes presence-absence atlas data and aggregates data to larger grain sizes, returning occupancy at each grain size for use in  \code{\link{downscale}} modelling. Atlas data may be in the form of raster data or as a data frame with presence-absence information and cell coordinates.

The extent for all scales is standardised to that of the largest grain size by applying a threshold for the proportion of unsampled atlas cells allowed within a cell at the largest grain size. The threshold can be chosen by the user or one of four threshold selections methods apply. See \code{\link{upgrain.threshold}} for a visualisation of the thresholds.

The function outputs a data frame of occupancies suitable as input for \code{\link{downscale}}, and will also plot the original atlas data along with the standardised data for each upgrained scale. In all plots presence = red, absence = white, and NA = grey.
}
\details{
When aggregating data to larger grain sizes, the extent of all grain sizes is set to the extent of the largest grain size. At the atlas scale, unsampled NA cells that fall within this extent are assigned as absences. This ensures that there are subsequently equal extents across all scales and therefore consistency of occupancy estimates. Therefore, it is necessary to apply a threshold whereby only those cells at the largest grain size are retained if a certain proportion of them are sampled at the altas scale. However, there is a trade-off between assigning unsampled cells as absences, and discarding sampled cells and presences.

\if{html}{It is highly advisable before selecting a threshold to explore this trade-off on a case-by-case basis with \code{\link{upgrain.threshold}} and to read the help file therein as well as referring to \code{vignette("Upgraining", package = "downscale")}.}

\if{latex}{It is highly advisable before selecting a threshold to explore this trade-off on a case-by-case basis with \code{\link{upgrain.threshold}} and to read the help file therein as well as referring to:\cr
\code{vignette("Upgraining", package = "downscale")}.}

The user may input their own threshold or use one of four suggested threshold criteria:
\if{html}{
\tabular{llll}{
      Threshold \tab  Method \tab \tab Description \cr
      0  \tab All_Sampled	\tab \tab All of the original atlas data is included. \cr
      Species specific \tab	All_Occurrences \tab \tab	The threshold where no occurrences in the atlas data are       excluded. \cr
      Atlas specific \tab	Gain_Equals_Loss \tab \tab	The threshold where the number of sampled atlas cells reclassified as No Data equals the number of unsampled exterior cells reclassified as absence. In this threshold the new standardised extent also equals the extent of the original atlas data. \cr
      1 \tab	Sampled_Only \tab \tab	Only cells that contain 100\% sampled atlas data are included. \cr
      }
}
\if{latex}{\figure{Threshold.png}{options: width=10cm}}

\if{html}{\figure{Maps.png}{options: width="75\%" alt="Figure: Maps.png"}}
\if{latex}{\figure{Maps.png}{options: width=13cm}}
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

## use a specified threshold - method must equal NULL
upgrain(atlas.data = atlas.data,
        cell.width= 10,
        scales = 3,
        threshold = 0.15,
        method = NULL)
        
## use one of the suggested methods - do not specify a threshold
all_sampled <- upgrain(atlas.data = atlas.data,
                       cell.width= 10,
                       scales = 3,
                       method = "All_Sampled")
        

## input data for downscale for All Sampled threshold
all_sampled$occupancy.stand

## and the original occupancies (note that extent varies with scale)
all_sampled$occupancy.orig

## if data is SpatialPointsDataFrame
presence <- data.frame(presence = atlas.data$presence)
spdf <- SpatialPointsDataFrame(coords = data.frame(x = atlas.data$x, 
                                                   y = atlas.data$y),
                               data = presence)
                               
all_sampled <- upgrain(atlas.data = spdf,
                       cell.width= 10,
                       scales = 3,
                       method = "All_Sampled")  
}
\author{
Charles Marsh <\email{charliem2003@gmail.com}>
}
\references{
Groom, Q., Marsh, C.J., Gavish, Y. Kunin, W.E. (2018). How to predict fine resolution occupancy from coarse occupancy data, \emph{Methods in Ecology and Evolution}. In press.

Marsh, C.J, Barwell, L.J., Gavish, Y., Kunin, W.E. (2018). downscale: An R package for downscaling species occupancy from coarse-grain data to predict occupancy at fine-grain sizes", \emph{Journal of Statistical Software, Code Snippets} 86(3), 1-20.
}

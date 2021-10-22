\name{upgrain.threshold}
\alias{upgrain.threshold}
\title{Exploration of trade-offs in threshold selection for upgraining}
\usage{
upgrain.threshold(atlas.data, cell.width = NULL, scales,
                  thresholds = seq(0, 1, 0.01))
}
\arguments{
\item{atlas.data}{either a raster file of presence-absence; an object of class SpatialPointsDataFrame atlas data where the data frame contains a single column of presence (1) and absences (0); or a data frame of sampled cells, which must contain columns named 'x', 'y' and 'presence'.}

\item{cell.width}{if data is a data frame or SpatialPointsDataFrame, the cell widths of sampled cells. If data is a raster then leave as default (\code{= NULL})}

\item{scales}{the number of cells to upgrain. Upgraining will happen by
  factors of 2 - ie if scales = 3, the atlas data will be aggregated in 2x2
  cells, 4x4 cells and 8x8 cells.}
  
\item{thresholds}{a vector of thresholds between and including 0 and 1 for the quantity of unsampled NA cells that can be included.}
}
\value{
Returns a list contatining two objects:
\item{Thresholds}{the threshold values for the four default threshold selections.}
\item{Data}{Data frame containing six columns:
   \tabular{lll}{
     \code{Thresholds} \tab  \tab Thresholds tested.\cr
     \code{SampledExluded} \tab  \tab Number of sampled cells exluded.\cr
     \code{SampledIncluded} \tab  \tab Number of sampled cells included.\cr
     \code{UnsampledAdded} \tab  \tab Number of unsampled NoData cells added.\cr
     \code{Extent} \tab  \tab Total number of cells included.\cr
     \code{OccurrencesExluded} \tab  \tab Number of cells with presence records excluded.\cr
     }
}
}
\description{
Explores the NoData threshold selection for upgraining whilst keeping a constant extent across scales. The thresholds are the quantity of unsampled cells at the atlas scale allowed within each cell at the largest grain size. A low threshold means that many unsampled cells will be assigned as absences, whereas a high threshold will mean that many sampled cells and many presence records will be excluded. These trade-offs are plotted, and four possible threshold choices are suggested and their maps presented.
}
\details{

A more detailed description is available at \code{vignette("Upgraining", package = "downscale")}.

In order to \code{\link{downscale}} we need to \code{\link{upgrain}} our atlas data across several scales. However, if the atlas data is not rectangular, as we aggregate cells during upgraining then the extent also increases.

\if{html}{\figure{Original.png}{options: width="100\%" alt="Figure: Original.png"}}
\if{latex}{\figure{Original.png}{options: width=14cm}}

Instead we must ensure the extent is constant across all scales by fixing the extent at all grain sizes to the extent of the largest grain size and convert our proportion of occupied cells back to area of occupancy by using the standardised extent (not the original atlas data extent).

However, if we fix the extent there is trade-off between assigning large areas of unsampled areas as absence, and discarding sampled areas and known Occurrences. The \code{upgrain.threshold} function allows visualisations of this trade-off at the atlas scale through four plots:

\if{html}{
\tabular{lll}{
a) \tab \tab  The total standardised extent; \cr
b) \tab \tab  The number of unsampled cells added and assigned as absences, and the number of sampled cells excluded and assigned as No Data; \cr
c) \tab \tab  The proportion of the original atlas data retained; \cr
d) \tab \tab  The proportion of known Occurrences excluded. \cr
}
}
\if{latex}{
\tabular{lll}{
a) \tab \tab  The total standardised extent; \cr
b) \tab \tab  The number of unsampled cells added and assigned as absences, \cr
   \tab \tab and the number of sampled cells excluded and assigned as No Data; \cr
c) \tab \tab  The proportion of the original atlas data retained; \cr
d) \tab \tab  The proportion of known Occurrences excluded. \cr
}
}
The final choice of threshold is up to the user on a case-by-case basis but we propose four threshold criteria in this function:
\if{html}{
\tabular{llll}{
      Threshold \tab  Method \tab \tab Description \cr
      0  \tab All_Sampled  \tab \tab All of the original atlas data is included. \cr
      Species specific \tab	All_Occurrences \tab \tab	The threshold where no occurrences in the atlas data are       excluded. \cr
      Atlas specific \tab	Gain_Equals_Loss \tab \tab	The threshold where the number of sampled atlas cells reclassified as No Data equals the number of unsampled exterior cells reclassified as absence. In this threshold the new standardised extent also equals the extent of the original atlas data. \cr
      1 \tab	Sampled_Only \tab \tab	Only cells that contain 100\% sampled atlas data are included. \cr
      }
}
\if{latex}{\figure{Threshold.png}{options: width=10cm}}

The function also creates maps for each of these four thresholds. In the example case this clearly demonstrates the trade-off between generating assumptions about unsampled areas, and losing data (and Occurrences) for the sampled atlas data.

\if{html}{\figure{Maps.png}{options: width="100\%" alt="Figure: Maps.png"}}
\if{latex}{\figure{Maps.png}{options: width=13cm}}
}

\examples{
## example species data
data.file <- system.file("extdata", "atlas_data.txt", package = "downscale")
atlas.data <- read.table(data.file, header = TRUE)

## if the input data is a data frame it must have the columns "x", "y"
## and "presence"
head(atlas.data)

thresh <- upgrain.threshold(atlas.data = atlas.data,
                            cell.width = 10,
                            scales = 3,
                            thresholds = seq(0, 1, 0.02))

## the four optional thresholds
thresh$Thresholds
head(thresh$Data)
}
\author{
Charles Marsh <\email{charliem2003@gmail.com}>
}
\references{
Groom, Q., Marsh, C.J., Gavish, Y. Kunin, W.E. (2018). How to predict fine resolution occupancy from coarse occupancy data, \emph{Methods in Ecology and Evolution}. In press.

Marsh, C.J, Barwell, L.J., Gavish, Y., Kunin, W.E. (2018). downscale: An R package for downscaling species occupancy from coarse-grain data to predict occupancy at fine-grain sizes", \emph{Journal of Statistical Software, Code Snippets} 86(3), 1-20.
}

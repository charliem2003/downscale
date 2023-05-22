# downscale: Downscaling Species Occupancy

An R package that provides a set of functions that model the occupancy-area relationship (OAR) of known coarse scale data. The models are then extrapolated to predict the proportion of occupied area at finer grain sizes.

To install the package directly from github (rather than CRAN) run:
`remotes::install_github("charliem2003/downscale", build_vignettes = TRUE)`

### Overview

The package provides three sets of functions for each stage of analysis:

1) `upgrain` and `upgrain.threshold` prepare atlas data for downscaling. 

2) `downscale` and `hui.downscale` model the OAR to the prepared data for one of ten possible downscaling models. 

3) `predict.downscale` and `plot.predict.downscale` take the model outputs and predict occupancy at finer grains. 

Finally, `ensemble.downscale` will run `downscale` and `predict.downscale` for a number of selected downscaling functions and calculate the mean predicted occupancies across all models.

The general flow of the package, and the inputs required for each function, is as follows:

<img src="https://user-images.githubusercontent.com/46057079/143235954-b7e1c851-142c-4e95-bcf5-485024a1f75d.png" width="500" height="600">

Two vignettes are available to guide users. Both work through examples in code:

`vignette("Downscaling", package = "downscale")`

`vignette("Upgraining", package = "downscale")`

Or are available through the github wiki:

[Introduction to downscaling species occupancy](https://github.com/charliem2003/downscale/wiki/Introduction-to-downscaling-species-occupancy)

[Upgraining atlas data for downscaling](https://github.com/charliem2003/downscale/wiki/Upgraining-atlas-data-for-downscaling)

### Credits

This package was created as part of deliverable D3.2 of WP3 of the project: **EU-BON: Building the European Biodiversity Observation Network} - a 7th Framework Programme funded by the European Union under Contract No. 308454.**

Author: Charles Marsh with input from Louise Barwell and Cang Hui.
Maintainer: Charles Marsh <charlie.marsh@mailbox.com>
Website: https://github.com/charliem2003/downscale

For reporting bugs or requesting information please include *'downscale'* in the subject line.

## References

Azaele, S., Cornell, S.J., & Kunin, W.E. (2012). Downscaling species occupancy from coarse spatial scales. *Ecological Applications* 22, 1004-1014.

Barwell, L.J., Azaele, S., Kunin, W.E., & Isaac, N.J.B. (2014). Can coarse-grain patterns in insect atlas data predict local occupancy? *Diversity and Distributions* 20, 895-907.

Hui, C. (2009). On the scaling patterns of species spatial distribution and association. *Journal of Theoretical Biology* 261, 481-487.

Hui, C., McGeoch, M.A., & Warren, M. (2006). A spatially explicit approach to estimating species occupancy and spatial correlation. *Journal of Animal Ecology* 7, 140-147.
  
Groom, Q., Marsh, C.J., Gavish, Y. Kunin, W.E. (2018). How to predict fine resolution occupancy from coarse occupancy data. *Methods in Ecology and Evolution* 9(11), 2273-2284.

Marsh, C.J, Barwell, L.J., Gavish, Y., Kunin, W.E. (2018). downscale: An R package for downscaling species occupancy from coarse-grain data to predict occupancy at fine-grain sizes. *Journal of Statistical Software* 86(Code Snippet 3), 1-20.

Marsh, C.J, Gavish, Y., Kunin, W.E., Brummitt N.A. (2019). Mind the gap: Can downscaling Area of Occupancy overcome sampling gaps when assessing IUCN Red List status? *Diversity and Distributions* 25, 1832-1845.

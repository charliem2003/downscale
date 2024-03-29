# downscale 5.0.0
* 24/05/2023
* SWITCHING TO terra + sf FROM raster + sp
* - Simple reformatting
* - optimisation functions renamed to optimisePars**
* - email address updated
* - minor changes to plots
* - error checking refined and moved to single function

# downscale 4.2-1
* 24/11/2021
* - Wesbite added to namespace

# downscale 4.2-0
* 24/11/2021
* downscale resubmission to CRAN. Amongst other minor changes:
* - vignettes reworked in markdown
* - new refs and a citation added
* - on.exit used for resetting par() in plotting function
* - exportClasses removed from namespace

# downscale 3.0-1
* 31/08/2018
* - Extra information and reference added to description

# downscale 3.0-0
* 30/08/2018
* - CITATION file added for JSS article and MEE article added to references

# downscale 2.0-3
* 24/08/2017
* - CITATION file added for JSS in press article

# downscale 2.0-2
* 08/08/2017
* - SpatialPointsDataFrame input added to upgrain helpfile
* - downscale helpfile updated

# downscale 2.0-1
* 01/08/2017
* - Bug when using SpatialPointsDataFrame corrected

# downscale 2.0-0
* 28/07/2017
* - 'lat' and 'lon' replaced with 'x' and 'y' in data frame input
* - SpatialPointsDataFrame allowed as input

# downscale 1.2-4
* 24/11/2016
* - Bug fixed with error checking in predict.downscale and ensemble.downscale

# downscale 1.2-3
* 29/09/2016
* - No need to input cell width in hui.downscale with raster or upgrain object

# downscale 1.2-2
* 18/08/2016
* - Default upgraining method changed to "All_Sampled"
* - Names of parameters harmonised with JSS paper

# downscale 1.2-1
* 10/03/2016
* - Tutorials checked and adjusted

# downscale 1.2-0
* 09/03/2016
* - 'spocc' replaced by 'rgbif' in tutorial

# downscale 1.1-0
* 25/02/2016
* - Bug for calculating extent in upgrain function fixed
* - upgrain now can return all upgrained rasters
* - ExtendRaster function deleted (now obsolete)

# downscale 1.0
* 09/10/2015
* - upgrain and upgrain.threshold: 'All Presences' changed to 'All_Occurrences'
* - upgrain and upgrain.threshold: bugs fixed when selecting number of scales
* - upgrain: plotting now optional
* - DataInput: bug in scale of endemism corrected
* - Vignettes updated
* - plot and predict changed to s3 methods
* - help files and examples updated
* - Tutorial vignette added
* - Raster package now depends
* - Argument names harmonised across functions
* - Option to specify starting parameters
* - More plotting arguments
* - Calculates AOO as well as occupancy
* - Upgraining vignette added

# downscale 0.7
* 24/04/2015
* - Function upgrain.threshold added to explore trade-offs in threshold selection
* - upgrain now allows thresholds of the four default thresholds

# downscale 0.6
* 14/04/2015
* - Plotting added to upgrain function
* - Hui model added

# downscale 0.5
* 13/03/2015
* - Starting parameters modified
* - Corrected plotting when 0's predicted

# downscale 0.5
* 24/02/2015
* - ensemble.predict: means calculated as mean of log occupancies
* - ensemble.predict: warning messages for inconsistent results
* - package.rd file updated
* - ensemble * - different tolerances allowed for modelling and predicting

# downscale 0.4
* 05/02/2015
* - plot=TRUE added to ensemble plotting
* - optimisation of logistic model now includes lower bounds
* - optimisation of GNB model now includes upper and lower bounds

# downscale 0.3
* 04/02/2015
* - ensemble.downscale function added

# downscale 0.2
* 03/02/2015
* - plotting function added and plotting option added to 'predict.downscale'
* - output of function downscale defined as class 'downscale'
* - output of function predict.downscale defined as class 'predict.downscale'
* - outputs updated to include observed data for use in the plot function
* - help files for downscale and predict.downscale updated
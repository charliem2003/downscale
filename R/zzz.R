.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
  cat("From v5.0.0 objects of class 'raster' and 'sp' are no longer supported.
  Please use 'terra' or 'sf' inputs instead.
  You can convert to 'terra' or 'sf' objects through rast() or st_as_sf()"))
}
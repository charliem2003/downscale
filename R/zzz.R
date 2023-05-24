.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "From v5.0.0 objects of class 'raster' and 'sp' are being phased out.
They are still supported but we recommend using 'terra' or 'sf' inputs instead.
You can convert to 'terra' or 'sf' objects through rast() or st_as_sf()"
  )
}
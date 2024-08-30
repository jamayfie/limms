#' Calculate parts per million (ppm)
#'
#' This function calculates m/z deviations in parts per millon (ppm).
#'
#' @param	expected	The expected m/z value.
#' @param	observed 	The measured m/z value.
#' @return	Returns the deviation from the expected m/z in parts per million.
#' @details	A simple ppm calculator and internal helper function in limms
#' @export 
#' @seealso \code{\link{dbMatch}}
#' @examples
#' # Is a [M+H]+ measured at 90.0552 alanine?  
#'
#' ppmCalc(90.0550, 90.0552)
#' # May not be alanine, but with a ppm < 2.23 its close.

ppmCalc <- function(expected, observed) {((observed-expected)/expected)*1000000}

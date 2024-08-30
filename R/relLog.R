#' Calculate relative log expression (RLE)
#'
#' This function calculates relative log expression (RLE) for QC plotting.
#' @param	x	An R object containing the data to be scaled arranged in
#' columns.
#' @param	intensities	The column numbers containing the intensities to be transformed.
#' The default is to assume all columns are intensities.
#' @return	Returns a data.table with relative log intensities appended to the input.
#' @details	This function transforms columns of data to relative log expression.
#' An internal function generally used in conjunction with boxplots to aid
#' the display of data that contain many measurements.
#' @note	Works on any numeric column, but data should already be log transformed.
#' @import data.table
#' @export
#' @seealso \code{\link{runNorm}}
#' @examples
#' 
#' require(data.table)
#' 
#' # A simplified version of the example from the runNorm function.
#'
#' # For the example CBS dataset CBS.xcms_diffreport,
#' # An xcms diffreport included as limms package data
#'
#' # Impute zeros and log transform:
#' all.log <- imputeZerosUnifMin(CBS.xcms_diffreport[,24:55])
#' 
#' # define colors for plots
#' colMet <- c("palegreen", "darkgreen")[c(rep(2,16), rep(1,8), rep(2,8))]
#'
#' # Boxplots of log intensities
#' par(mfcol=c(1,2), oma=c(4,1,1,1))
#' boxplot(all.log[, .SD, .SDcols=patterns("log2")], col=colMet, las=2,  
#'   main="Unormalized intensities (log)")
#'
#' # Boxpots of Relative log expression (RLE)
#' CBS.relLog <- relLog(all.log, intensities=33:64)
#' 
#' boxplot(CBS.relLog[, .SD, .SDcols=patterns("rle")], col=colMet, las=2, ylim=c(-5, 5),
#'   main="RLE unormalized intensities")
#' abline(h=0)

relLog <-
  function(x, intensities = 1:nrow(x))
  {
    ms.events <- . <- NULL
    if (is.data.table(x) == TRUE) {
      y <- copy(x)
    }
    else {
      y <- data.table(x)
    }
    
    event.names <- "ms.events"
    y[, ms.events := 1:nrow(y)]
    int.names <- names(y[, ..intensities])
    
    yrle <-
      y[, .(ms.events, 
        t(scale(t(.SD), center = unlist(y[, .(apply(.SD, 1, stats::median)), .SDcols = int.names]), 
        scale = FALSE))), .SDcols = int.names]
    setnames(yrle, int.names, paste0(int.names, ".rle"))

    y[yrle, on = event.names][, ms.events := NULL]
    
    
  }



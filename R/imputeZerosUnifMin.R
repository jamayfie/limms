#' Imputation of zeros using random minima
#'
#' This function addreses zero values by imputation.
#' For convenience, it can also log2 transform the data.
#' @param	x	An R object containing the data to be imputed arranged in
#' columns.
#' @param	output	Either "log2" or "imputed".  Defaults to "log2".
#' @param	intensities	Column range containing the mass spectrometry peak intensities. 
#' Enter the columns containing intensity data by number, for example 24:69 in the CBS data,
#' or enter by column names.  
#' Defaults to using all columns, which only works if the input contains only intensity data.
#' @param	seed	Change the seed for consistent output.
#' @return	If the output is set to "log2", the default, x is returned
#' as a data table plus appended columns with a ".log2" suffix for each
#' column specified in the "intensities" argument.
#' Zeros are replaced and intensity values log2 transformed in these columns.
#' If the output is set to "imputed", x is returned with columns appended 
#' and a ".imputed" suffix added for each column in the "intensities" argument.
#' Zeros are replaced in these columns without transformation.
#' @details	The minimum value within each column is calculated.
#' Zero values are replaced at random by runif, such that
#' imputed values are between 0 and the column minimum.
#' Log transformation of data is preferred for downstream analysis, and is
#' implemented here for convenience.
#' Minimum intensities < 1 are possible; hence, to avoid negative log2 values,
#' 1 is added to all intensities prior to log2 transformation.
#' @note	Due to the random replacement of zeros, imputed values will change 
#' each time the function is called.  The seed argument calls set.seed() 
#' to maintain consistent values, but should be changed from the default.
#' @seealso \code{\link{runif}}
#' \code{\link{runNorm}}
#' @family limms
#' @import data.table
#' @export 
#' @examples
#' 
#' require(data.table)
#' 
#' # For the CBS dataset included as limms package data, 
#' # the xcms diffreport CBS.xcms_diffreport, 
#' # choose columns with data, for example, 24 to 55 contain non-spiked samples.
#' 
#' # Return imputed values using:
#' imputeZerosUnifMin(CBS.xcms_diffreport, intensities=24:55, seed=890, output="log2")
#'
#' # Return lists of imputed values and minima.
#' imputeZerosUnifMin(CBS.xcms_diffreport, intensities=24:55, seed=890, output="impounded")
#'
#' # To store imputed, log transformed values, call to a new object
#' # Change the seed, 478, to a different number for each new data set.
#' all.i_l <- imputeZerosUnifMin(CBS.xcms_diffreport, intensities=24:55, seed=478)
#'
#' # The imputed measurements can be flagged and counted.
#' metimp <- t(sapply(1:dim(CBS.xcms_diffreport)[1], function(i) CBS.xcms_diffreport[i, 24:55]==0))
#' numimp <- sapply(1:dim(CBS.xcms_diffreport)[1], function(i) length(which(metimp[i,]=="TRUE")))
#' 
#' # Quick QC according to the number of impounded measurements
#' # which peaks have no samples with impounded zeros?
#' summary(numimp == 0)
#' # which peaks have > 5 impounded?
#' summary(numimp > 5)
#'
#' lcids <- sapply(1:dim(CBS.xcms_diffreport)[1], function(i) 
#'   paste(as.character(
#'     colnames(CBS.xcms_diffreport[, 24:55])[which(CBS.xcms_diffreport[i, 24:55]==0)]),
#'     collapse="; "))
#' 
#' # Generate a table with information about impounded measurements
#' # including how many and which samples were impounded for each metabolite
#' xi <- data.table(cbind(numimp,lcids))
#' head(xi)



imputeZerosUnifMin <-
  function(x,
           output = "log2",
           intensities = 1:ncol(x),
           seed = 99)
  {
    ms.events <- . <- NULL
    set.seed(seed)
    if (is.data.table(x) == TRUE) {
      xc <- copy(x)
    }
    else {
      xc <- data.table(x)
    }
    event.names <- "ms.events"
    xc[, ms.events := 1:nrow(xc)]
    int.names <- names(xc[, ..intensities])
    
    m <- copy(xc[, ..intensities])
    for (j in int.names)
      set(m, j = j, value = min(m[[j]][m[[j]] != 0]))
    m <- as.numeric(unlist(m[1]))
    
    xm <- xc[, .SD, .SDcols = c(event.names, int.names)]
    for (j in int.names)
      set(xm, which(xm[[j]] == 0), j, stats::runif(sum(xm[[j]] == 0), 0, m))
    if (output == "log2")	 {
      for (j in int.names)
        set(xm, j = j, value = log2(xm[[j]] + 1))
      setnames(xm, old = int.names, new = paste0(int.names, ".log2"))
    }
    else if (output == "imputed") {
      setnames(xm, old = int.names,
               new = paste0(int.names, ".imputed"))
    }
    
    xc[xm, on = event.names][, ms.events := NULL]
  }	
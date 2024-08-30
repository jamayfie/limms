#' Normalization of sample data
#'
#' This function provides methods for within sample normalization of data.
#' The suggested method is full quantile normalization, but global normalization
#' or principal component normalization are also available.
#' @param	x	An R object containing data to be normalized.
#' @param	intensities	Numeric column range.  Enter the columns containing intensity
#' data by number, for example 24:69 in the CBS data set, or enter by column names.
#' Using column names via pattern match to the "log2" suffix can be an advantage 
#' in the future or when re-using code.
#' Defaults to using all columns, which only works if the input contains only intensity data.
#' @param	method	One of 4 normalization methods, 
#' "fullQuantileNorm" (default), "globalScalingNorm", "pc1Norm", or a 
#' user defined function may be used by specifying "custom".
#' @param	FUN	If method is set to "custom", the FUN argument names the
#' function used to normalize the data
#' @details	The object x must be arranged with each sample as a column.
#' Zero values are not tolerated by "fullQuantileNorm", and should be removed
#' prior to normalization using imputeZerosUnifMin.
#' @seealso \code{\link{rle}}
#' @family limms
#' @import aroma.light
#' @import data.table
#' @import ggplot2
#' @import forcats
#' @import reshape2
#' @export
#' @examples
#'
#' require(data.table)
#' require(aroma.light)
#' require(ggplot2)
#' require(forcats)
#' require(reshape2)
#'
#' # Select the data columns in the CBS data set, CBS.xcms_diffreport
#' # impute the zeros, and log2 transform (the default in imputeZerosUnifMin).
#' all.i_l <- imputeZerosUnifMin(CBS.xcms_diffreport, intensities=24:55, seed=478)
#'
#' # Normalization by full quantiles
#' all.fq <- runNorm(all.i_l, intensities=70:101)
#'
#' # Normalization by global scaling
#' all.gs <- runNorm(all.i_l, intensities=70:101, method="globalScalingNorm")
#'
#' # Normalization by principal components
#' all.pc1 <- runNorm(all.i_l, intensities=70:101, method="pc1Norm")
#'
#' # Normalization using a custom function, the sample median
#' # For demonstration only
#' medPeak <- function(x) {(x*apply(x,2,median))/x}
#' all.med <- runNorm(all.i_l, intensities=70:101, method="custom", FUN=medPeak)
#'
#'
#' # What did the normalization do?
#'
#' # add a color key
#' desMetB6 <- cbind(names(CBS.xcms_diffreport[,24:55]), c(rep("CBS",4), rep("CBS",4), rep("CBS",4), 
#' rep("CBS",4), rep("CBS",4), rep("CBS",4), rep("G307S",4), rep("G307S",4)), c(rep("Yes",4), 
#' rep("Yes",4), rep("Yes",4), rep("Yes",4), rep("No",4), rep("No",4), rep("Yes",4), 
#' rep("Yes",4)),  c(rep("High",4), rep("High",4), rep("Low",4), rep("Low",4), rep("High",4), 
#' rep("Low",4), rep("High",4), rep("Low",4)), c(rep(2,4), rep(1,4), rep(2,4), rep(1,4), 
#' rep(1,4), rep(1,4), rep(2,4), rep(2,4)))
#' desMetB6 <- data.frame(desMetB6)
#' names(desMetB6) <- c("Run", "Strain", "Met", "B6", "Rep")
#'
#' # methionine starvation induced the largest changes
#' # make the methionine starvation condition pale green
#' colCBS <- data.table(desMetB6)
#'
#'
#' # Boxplots of log intensities
#'
#' all.fq.melt <- data.table(reshape2::melt(all.fq, id.vars=1:69, value.name="intensity", 
#'  variable.name="sample"))
#' all.fq.melt[, Run := gsub("\\..*$", "", sample)]
#' all.fq.melt[, transform := gsub("^.*\\.", "", sample)]
#' all.fq.melt[, transform := gsub("norm", "quantile", transform)]
#'
#' all.fq.melt <- merge(all.fq.melt, colCBS, by="Run", sort=FALSE)
#'
#' ggplot(all.fq.melt, aes(x=fct_inorder(Run), y=intensity)) +
#'   geom_boxplot(aes(fill=Met)) +
#'   facet_wrap(~transform, scales="free_x") +
#'   theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
#'   xlab("sample") +
#'   ggtitle("comparison of log2 transformed to full quantile normalized data")
#'
#'
#' # Comparison of normalized data
#' norm.list <- list(quantile=all.fq[, 102:133], global_scaling=all.gs[, 102:133], 
#'  pc1=all.pc1[, 102:133], median=all.med[, 102:133])
#' all.norm.melt <- data.table(reshape2::melt(norm.list, value.name="intensity", 
#'  variable.name="sample"))
#'
#' all.norm.melt[, Run := gsub("\\..*$", "", sample)]
#' all.norm.melt[, transform := gsub("^.*\\.", "", sample)]
#' all.norm.melt[, transform := gsub("norm", "quantile", transform)]
#' all.norm.melt <- merge(all.norm.melt, colCBS, by="Run", sort=FALSE)
#'
#' ggplot(all.norm.melt, aes(x=fct_inorder(Run), y=intensity)) +
#'   geom_boxplot(aes(fill=Met)) +
#'   facet_wrap(~L1, scales="free") +
#'   theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
#'   xlab("sample") +
#'   ggtitle("comparison of normalization methods")
#'
#'
#' # Boxpots of Relative log expression (RLE)
#' # Normalizations like PC1 have different y-axis scales, making comparisons to
#' # other methods difficult. relLog is a limms internal function useful for these cases.
#'
#' norm.list <- list(quantile=relLog(all.fq, intensities=102:133)[, .SD, .SDcols=patterns("rle")],
#' global_scaling=relLog(all.gs, intensities=102:133)[, .SD, .SDcols=patterns("rle")],
#' pc1=relLog(all.pc1, intensities=102:133)[, .SD, .SDcols=patterns("rle")],
#' median=relLog(all.med, intensities=102:133)[, .SD, .SDcols=patterns("rle")])
#' all.norm.melt <- data.table(reshape2::melt(norm.list, value.name="intensity", 
#'  variable.name="sample"))
#' all.norm.melt[, Run := gsub("\\..*$", "", sample)]
#' all.norm.melt <- merge(all.norm.melt, colCBS, by="Run", sort=FALSE)
#'
#' ggplot(all.norm.melt, aes(x=fct_inorder(Run), y=intensity)) +
#'   geom_boxplot(aes(fill=Met)) +
#'   facet_wrap(~L1, scales="free") +
#'   theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
#'   xlab("sample") +
#'   ggtitle("comparison of normalization by relative log expression")
#'


runNorm <-
  function(x,
           method = "fullQuantileNorm",
           intensities = 1:ncol(x),
           FUN = NA)
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
    ynorm <- y[, ..intensities]
    
    if (method == "fullQuantileNorm") {
      . <- NULL
      fullQuantileNorm <- function(z, int = int.names)
      {
        z[, .(aroma.light::normalizeQuantile(as.matrix(z[, int, with = FALSE]), robust = TRUE))]
      }
      ynorm <- fullQuantileNorm(y)
    }
    
    else if (method == "globalScalingNorm") {
      . <- NULL
      globalScalingNorm <- function(z)
      {
        scale(z, center = (unlist(z[, .(lapply(.SD, function(i)
          stats::quantile(unlist(i, use.names = FALSE), probs = 0.75))), .SDcols = names(z)]) - 
            unlist(z[, .(stats::quantile(unlist(.SD, use.names = FALSE), probs = 0.75)),
            .SDcols = names(z)], use.names = FALSE)), scale = FALSE)
      }
      ynorm <- y[, .(globalScalingNorm(ynorm))]
    }
    
    else if (method == "pc1Norm") {
      . <- NULL
      pc1Norm <- function(z)
      {
        cz <- copy(z)
        s <- cz[, .(svd(.SD)$u)]
        s <- s[, 1]
        for (j in int.names) {
          set(cz,
              j = j,
              value = cz[[j]] - sum(cz[[j]] * s) * s / sum(s ^ 2))
        }
        cz
      }
      ynorm <- pc1Norm(ynorm)
    }
    
    else if (method == "custom") {
      ynorm <- FUN(ynorm)
    }
    
    setnames(ynorm, names(ynorm), paste0(names(ynorm), ".norm"))
    ynorm[, (event.names) := y[, .SD, .SDcols = event.names]]
    y[ynorm, on = event.names][, ms.events := NULL]
  }

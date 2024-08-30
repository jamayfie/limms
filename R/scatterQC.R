#' Plot mass spectrometry data for quality control
#'
#' This function makes two plots from a table of mass spectrometry data to assess quality
#' and to compare sample types.  Data are log transformed: input data should be raw,
#' untransformed peak data. 
#' 
#' 
#' The top scatter plot of nonzero ions shows the number of peaks measured above threshold
#' for each sample, binned by sample class.  
#' Nonzero ions are counted first to make the top plot, but then are imputed 
#' using 'imputeZerosUnifMin()' to generate the lower plot.
#' 
#' The bottom boxplot shows the quantiles of the imputed peak intensities.
#' 
#' These plots are useful compared to historical data and to compare the sample classes.
#' Note: if 'xcms' was used to generate the peaklistCols, 
#' 'fillPeaks()' adds data to zero peaks: avoid using 'fillPeaks()' to get a true nonzero
#' peak list, but a fair number of zeros will still remain after 'fillPeaks()'.
#'
#' @param	peaklistCols	An R object of peak data, with peaks in rows and samples in columns.
#' @param	seed 	Set a seed to generate the same plot each time the function is used.
#' @param groupNames	If used, a vector of names for each column with peak data.
#' To group the samples, enter a vector of group names as long as the column names.
#' If not used, the names attribute of the input object is used and no groups are plotted.
#' @param main	The main plot title.
#' @param outlayer	If a sample is an outlier according to 'check_outlier()', 
#' setting outlayer=TRUE adds a text label, allowing identification.  
#' The default is FALSE.
#' @param zero	Set to zero=TRUE to start the y axis at zero.
#' The default is FALSE, zooming the y-axis to the data.  
#' @import data.table
#' @import ggplot2
#' @import forcats
#' @import cowplot
#' @importFrom stats median
#' @importFrom grDevices boxplot.stats
#' @export 
#' @seealso \code{\link{dbMatch}}
#' \code{\link{imputeZerosUnifMin}}
#' @examples
#'
#' require(data.table)
#' require(ggplot2)
#' require(cowplot)
#' require(forcats)
#' 
#' # Extract the data columns from CBS.xcms_diffreport, 
#' # an xcms diffreport object included as limms package data,
#' # impute the zeros, and log2 transform (the default in imputeZerosUnifMin).
#'
#' desMetB6 <- cbind(names(CBS.xcms_diffreport[,24:55]), c(rep("CBS",4), rep("CBS",4), rep("CBS",4), 
#' rep("CBS",4), rep("CBS",4), rep("CBS",4), rep("G307S",4), rep("G307S",4)), 
#' c(rep("Yes",4), rep("Yes",4), rep("Yes",4), rep("Yes",4), rep("No",4), rep("No",4),
#' rep("Yes",4), rep("Yes",4)),  c(rep("High",4), rep("High",4), rep("Low",4), 
#' rep("Low",4), rep("High",4), rep("Low",4), rep("High",4), rep("Low",4)), c(rep(2,4),
#' rep(1,4), rep(2,4), rep(1,4), rep(1,4), rep(1,4), rep(2,4), rep(2,4)))
#' desMetB6 <- data.frame(desMetB6)
#' names(desMetB6) <- c("Run", "Strain", "Met", "B6", "Rep")
#'
#' # QC plot of all samples
#' # this can be useful for checking run order, etc
#' scatterQC(CBS.xcms_diffreport[, c(24:55)], seed=75, main="QC of all samples")
#'
#' # QC plot by group
#' # Treatments are expected to differ... but it's still useful to know which ones!
#' # flag outlier samples
#' scatterQC(CBS.xcms_diffreport[, c(24:55)], seed=75, groupNames=desMetB6$Met, 
#'   main="QC by methionine addition", outlayer=TRUE)
#' 


scatterQC <-
  function(peaklistCols,
           seed = NULL,
           groupNames = NULL,
           main = NULL,
           outlayer = FALSE,
           zero = FALSE) {
    outlier_Ions <- Ions <- Covariates <- label_Ions <- Sample <- . <- NULL
    Class <- Intensity <- Median <- median <- NULL
    patterns <- ggplot_buildggplot_build <- ggplot_build <- NULL
    Nions <- apply(peaklistCols, 2, function(i)
      sum(i != 0))
    
    if (is.null(seed)) {
      seed <- sample(1:1000, 1)
    }
    set.seed(seed)
    print("Warning: scatterQC performs a log transformation: input data should not be log transformed.")
    plCi <-
      imputeZerosUnifMin(peaklistCols)[, .SD, .SDcols = patterns("log2")]
    
    # quantiles of the intensity after imputation
    q5 <- apply(plCi, 2, function(i) {
      stats::quantile(i, probs = 0.05)
    })
    q95 <- apply(plCi, 2, function(i) {
      stats::quantile(i, probs = 0.95)
    })
    q25 <- apply(plCi, 2, function(i) {
      stats::quantile(i, probs = 0.25)
    })
    med <- apply(plCi, 2, median)
    q75 <- apply(plCi, 2, function(i) {
      stats::quantile(i, probs = 0.75)
    })
    
    
    # Build an object of data to plot
    ionsSam <-
      data.frame(
        Sample = names(Nions),
        Ions = Nions,
        q5,
        q25,
        Median = med,
        q75,
        q95
      )
    
    
    # Groupings
    if (is.null(groupNames)) {
      groupNames <- ionsSam$Sample
    }
    
    
    ionsSamG <- data.table(Covariates = groupNames, ionsSam)
    
    
    # Flag outliers for labeling samples in QC scatterplots
    # http://stackoverflow.com/questions/33524669/labeling-outliers-of-boxplots-in-r
    check_outlier <- function(v, coef = 1.5) {
      quantiles <- stats::quantile(v, probs = c(0.25, 0.75))
      IQR <- quantiles[2] - quantiles[1]
      res <- v < (quantiles[1] - coef * IQR) |
        v > (quantiles[2] + coef * IQR)
      return(res)
    }
    
    ionsSamG[, outlier_Ions := check_outlier(Ions), by = Covariates]
    ionsSamG[, label_Ions := ifelse(outlier_Ions, as.character(ionsSamG[, Sample]), "")]
    
    iOut <- ionsSamG[outlier_Ions == TRUE]
    
    # Plot ions
    plotIons <-
      ggplot2::ggplot(ionsSamG, aes(x = forcats::fct_inorder(Covariates), y = Ions)) +
      ggplot2::geom_hline(
        yintercept = stats::median(Nions),
        color = "grey",
        linewidth = 0.25
      ) +
      ggplot2::geom_hline(
        yintercept = stats::quantile(Nions, probs = 0.95),
        linetype = "dotted",
        color = "grey",
        linewidth = 0.25
      ) +
      ggplot2::geom_hline(
        yintercept = stats::quantile(Nions, probs = 0.05),
        linetype = "dotted",
        color = "grey",
        linewidth = 0.25
      ) +
      ggplot2::geom_hline(
        yintercept = stats::quantile(Nions, probs = 0.75),
        linetype = "dashed",
        color = "grey",
        linewidth = 0.25
      ) +
      ggplot2::geom_hline(
        yintercept = stats::quantile(Nions, probs = 0.25),
        linetype = "dashed",
        color = "grey",
        linewidth = 0.25
      ) +
      
      ggplot2::geom_jitter(
        aes(color = Covariates),
        width = 0.25,
        height = 0.1,
        size = 3
      ) +
      
      ggplot2::ggtitle(paste0(main, "\nNonzero ions")) +
      ggplot2::theme(
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none"
      ) +
      ggplot2::theme(plot.margin = unit(c(0.1, 0.1, 0.3, 0.1), "cm"))
    
    
    if (outlayer == TRUE) {
      plotIons <-
        plotIons + ggplot2::annotate(
          "text",
          x = iOut[, Covariates],
          y = iOut[, Ions],
          label = iOut[, label_Ions],
          hjust = -0.2
        )
    }
    
    if (zero == TRUE) {
      plotIons <- plotIons + ggplot2::expand_limits(y = 0)
    }
    
    
    # get the intensity quantiles
    pl.melt <-
      data.table::melt(
        plCi,
        measure.vars = patterns(".log2"),
        variable.name = "samples",
        value.name = "intensity"
      )
    
    classNames <- data.table(samples = names(plCi), class = groupNames)
    pl.melt <- pl.melt[classNames, on = "samples"]
    
    plQuantiles <-
      data.frame(
        Sample = as.character(pl.melt$samples),
        Class = as.character(pl.melt$class),
        Intensity = pl.melt$intensity
      )
    ylim1 <- boxplot.stats(plQuantiles$Intensity)$stats[c(1, 5)]
    
    
    # Plot quantiles
    plotQuantiles <-
      ggplot2::ggplot(plQuantiles, aes(
        x = forcats::fct_inorder(Class),
        y = Intensity,
        fill = Class
      )) +
      ggplot2::geom_hline(
        yintercept = mean(ionsSamG[, Median]),
        color = "grey",
        linewidth = 0.25
      ) +
      ggplot2::geom_hline(
        yintercept = mean(ionsSamG[, q5]),
        linetype = "dotted",
        color = "grey",
        linewidth = 0.25
      ) +
      ggplot2::geom_hline(
        yintercept = mean(ionsSamG[, q95]),
        linetype = "dotted",
        color = "grey",
        linewidth = 0.25
      ) +
      ggplot2::geom_hline(
        yintercept = mean(ionsSamG[, q25]),
        linetype = "dashed",
        color = "grey",
        linewidth = 0.25
      ) +
      ggplot2::geom_hline(
        yintercept = mean(ionsSamG[, q75]),
        linetype = "dashed",
        color = "grey",
        linewidth = 0.25
      ) +
      
      ggplot2::geom_boxplot(notch = TRUE, outlier.shape = NA) +
      ggplot2::coord_cartesian(ylim = ylim1 * c(0.95, 1.05)) +
      
      ggplot2::ggtitle("Intensity boxplots") +
      ggplot2::theme(
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor.x = element_blank(),
        legend.position = "none"
      ) +
      ggplot2::theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
    
    
    # put the plots together in a grid
    
    p1c <- ggplot_gtable(ggplot_build(plotIons))
    p2c <- ggplot_gtable(ggplot_build(plotQuantiles))
    
    maxWidth <- grid::unit.pmax(p1c$widths[2:3], p2c$widths[2:3])
    
    p1c$widths[2:3] <- maxWidth
    p2c$widths[2:3] <- maxWidth
    
    cowplot::plot_grid(
      p1c,
      p2c,
      align = "v",
      rel_heights = c(1, 2),
      ncol = 1
    )
    
  }





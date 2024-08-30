#' 'limma' topTable object of CBS metabolomic data
#' 
#' The metabolite data in the data set 'CBS.xcms_diffreport'
#' analyzed using the limms pipeline to generate a table of
#' adjusted P values based on an F-score after contrasts of
#' all different sample classes.
#' The dataset is the output of the 'limma' function 'topTable',
#' available in the 'limms' function 'limmaTest'
#'
#' @name ttF.CBS
#' @docType data
#' @author Jacob Mayfield
#' @format A data frame with 52 rows and 9 columns
#' @references \url{https://pubmed.ncbi.nlm.nih.gov/22267502/}
#' @keywords data
#'
#' @seealso \code{\link{CBS.xcms_diffreport}}
#' \code{\link{limmaTest}}
#' @examples
#' data(ttF.CBS)
"ttF.CBS"

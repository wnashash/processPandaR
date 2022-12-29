#' @title Filter Gene Expression File
#' @description A simple function for additional gene expression filtering
#' @param exp (Object) Name of expression data frame object generated using read_gene()
#' @param filt (Character) Filter option for additional filtering - either sum, var, or both
#' @param thresh Numeric value between 0 and 1 for variance filter with default 0.3
#' @return Returns filtered gene expression data frame to global env and writes to working directory
#' @author Walid Nashashibi (\url{https://github.com/wnashash/})
#' @examples
#' \dontrun{
#' library(processNetZoo)
#' expression <- read_gene('expression_test.csv','gene')
#' generate_histogram(expression,'both')
#' filtered <- filter_gene(expression,'sum')
#' }
#' @import genefilter
#' @import utils
#' @export
#'
filter_gene <- function(exp,filt,thresh=NULL) {

  if(filt=='sum') {

    expData <- exp[rowSums(exp) >= 10,]

  } else if(filt=='var') {

    if(is.null(thresh)) {
      thr <- 0.3
    } else {
      thr <- thresh
    }

    expData <- exp[genefilter::rowSds(as.matrix(exp)) > thr*rowMeans(as.matrix(exp)),]

  } else if(filt=='both') {

    expData <- exp[rowSums(exp) >= 10,]

    if(is.null(thresh)) {
      thr <- 0.3
    } else {
      thr <- thresh
    }

    expData <- expData[genefilter::rowSds(as.matrix(expData)) > thr*rowMeans(as.matrix(expData)),]

  } else {

    stop("Filter not supported. Please, choose between sum, var, or both.", call. = FALSE)

  }

  utils::write.table(expData, file="expression.txt",
                     row.names=TRUE, col.names=FALSE, sep="\t", quote=F)

  return(expData)

}


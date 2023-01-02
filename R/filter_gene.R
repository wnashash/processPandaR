#' @title Filter Gene Expression File
#' @description A simple function for additional gene expression filtering
#' @param geneExp (Object) Name of expression data frame object generated using read_gene()
#' @param filter (Character) Filter option for additional filtering - either sum, var, or both
#' @param threshold Numeric value between 0 and 1 for variance filter with default 0.1
#' @return Returns filtered gene expression data frame to global env and writes to working directory
#' @author Walid Nashashibi (\url{https://github.com/wnashash/})
#' @examples
#' \dontrun{
#' library(processNetZoo)
#' expression <- read_gene('extdata/expression_test.csv','gene')
#' generate_histogram(expression,'both')
#' filtered <- filter_gene(expression,'sum')
#' }
#' @import genefilter
#' @import utils
#' @export
#'
filter_gene <- function(geneExp,filter,threshold=NULL) {

  if(filter=='sum') {

    expData <- geneExp[rowSums(geneExp) >= 10,]

  } else if(filter=='var') {

    if(is.null(threshold)) {
      thr <- 0.1
    } else {
      thr <- threshold
    }

    expData <- geneExp[genefilter::rowSds(as.matrix(geneExp)) > thr*rowMeans(as.matrix(geneExp)),]

  } else if(filter=='both') {

    expData <- geneExp[rowSums(geneExp) >= 10,]

    if(is.null(threshold)) {
      thr <- 0.1
    } else {
      thr <- threshold
    }

    expData <- expData[genefilter::rowSds(as.matrix(expData)) > thr*rowMeans(as.matrix(expData)),]

  } else {

    stop("Filter not supported. Please, choose between sum, var, or both.", call. = FALSE)

  }

  utils::write.table(expData, file="expression.txt",
                     row.names=TRUE, col.names=FALSE, sep="\t", quote=F)

  return(expData)

}


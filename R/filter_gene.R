#' @title Filter Gene Expression File
#' @description A simple function for additional gene expression filtering
#' @param expr (Object) Name of expression data frame object generated using read_gene()
#' @param filter (Character) Filter option for additional filtering - either sum, var, or both
#' @param threshold Numeric value between 0 and 1 for variance filter with default 0.1
#' @param expr_file_out (Character) Filename of the base group to write out results to.  No file will
#' be written when NULL.
#' @return Returns filtered gene expression data frame to global env and writes to working directory
#' @author Walid Nashashibi (\url{https://github.com/wnashash/})
#' @examples
#' \dontrun{
#' library(processNetZoo)
#' exp_path <- system.file("extdata",
#'      "expression_test.csv",
#'      package = "processNetZoo",
#'      mustWork = TRUE)
#'
#' expression <- read_gene(exp_path,'gene')
#' generate_histogram(expression,'both')
#' filtered <- filter_gene(expression,'sum')
#' }
#' @import genefilter
#' @import utils
#' @export
#'
filter_gene <- function(expr,
                        filter,
                        threshold=NULL,
                        expr_file_out = 'expression.txt') {

  if(filter=='sum') {

    expData <- expr[rowSums(expr) >= 10,]

  } else if(filter=='var') {

    if(is.null(threshold)) {
      thr <- 0.1
    } else {
      thr <- threshold
    }

    expData <- expr[genefilter::rowSds(as.matrix(expr)) > thr*rowMeans(as.matrix(expr)),]

  } else if(filter=='both') {

    expData <- expr[rowSums(expr) >= 10,]

    if(is.null(threshold)) {
      thr <- 0.1
    } else {
      thr <- threshold
    }

    expData <- expData[genefilter::rowSds(as.matrix(expData)) > thr*rowMeans(as.matrix(expData)),]

  } else {

    stop("Filter not supported. Please, choose between sum, var, or both.", call. = FALSE)

  }

  if(!is.null(expr_file_out)){
    utils::write.table(expData, file = expr_file_out,
                       row.names=TRUE, col.names=FALSE, sep="\t", quote=F)
  }

  return(expData)

}

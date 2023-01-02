#' @title Generate Histogram for Variance Visualization
#' @description A simple function to generate histogram of row sums and row variance
#' @param geneExp (Object) Name of expression data frame object generated using read_gene()
#' @param filter (Character) Filter option for visualization - either sum, var, or both
#' @return Returns histogram plot of row sums or row variance or both
#' @author Walid Nashashibi (\url{https://github.com/wnashash/})
#' @examples
#' \dontrun{
#' library(processNetZoo)
#' expression <- read_gene('extdata/expression_test.csv','gene')
#' generate_histogram(expression,'sum')
#' }
#' @import grDevices
#' @import graphics
#' @import stats
#' @export
#'
generate_histogram <- function(geneExp,filter) {

  if(filter=='sum') {

    graphics::hist(rowSums(geneExp),
                   col=grDevices::rgb(1,0,0,0.5),
                   xlab="Gene Expression Sums", main="Gene Exp Sums Across Samples")

  } else if(filter=='var') {

    graphics::hist(apply(geneExp,1,stats::var),
                   col=grDevices::rgb(0,0,1,0.5),
                   xlab="Gene Expression Variance", main="Gene Exp Variance Across Samples")

  } else if(filter=='both') {

    graphics::par(mfrow=c(1,2), mar=c(4,4,1,0))

    graphics::hist(rowSums(geneExp),
                   col=grDevices::rgb(1,0,0,0.5),
                   xlab="Gene Expression Sums", main="Sum")

    graphics::hist(apply(geneExp,1,stats::var),
                   col=grDevices::rgb(0,0,1,0.5),
                   xlab="Gene Expression Variance", main="Variance")

  } else {

    stop("Filter not supported. Please, choose between sum, var, or both.", call. = FALSE)

  }

}


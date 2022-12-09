#' @title Generate Histogram for Variance Visualization
#' @description A simple function to generate histogram of row sums and row variance
#' @param exp (Object) Name of expression data frame object generated using read_gene()
#' @param filt (Character) Filter option for visualization - either sum, var, or both
#' @return Returns histogram plot of row sums or row variance
#' @author Walid Nashashibi (\url{https://github.com/wnashash/})
#' @examples
#' \dontrun{
#' library(processPandaR)
#' expression <- read_gene('expression_test.csv','gene')
#' gen_hist(expression,'both')
#' }
#' @import grDevices
#' @import graphics
#' @import stats
#' @export
#'
gen_hist <- function(exp,filt) {

  if(filt=='sum') {

    graphics::hist(rowSums(exp),
                   col=grDevices::rgb(1,0,0,0.5),
                   xlab="Gene Expression Sums", main="Gene Exp Sums Across Samples")

  } else if(filt=='var') {

    graphics::hist(apply(exp,1,stats::var),
                   col=grDevices::rgb(0,0,1,0.5),
                   xlab="Gene Expression Variance", main="Gene Exp Variance Across Samples")

  } else if(filt=='both') {

    graphics::par(mfrow=c(1,2), mar=c(4,4,1,0))

    graphics::hist(rowSums(exp),
                   col=grDevices::rgb(1,0,0,0.5),
                   xlab="Gene Expression Sums", main="Sum")

    graphics::hist(apply(exp,1,stats::var),
                   col=grDevices::rgb(0,0,1,0.5),
                   xlab="Gene Expression Variance", main="Variance")

  } else {

    stop("Filter not supported. Please, choose between sum, var, or both.", call. = FALSE)

  }

}


#' @title Normalize Expression File
#' @description A simple function to normalize expression data and export
#' @param expData (Object) Name of expression data frame object generated using read_clean()
#' @param method (Character) Either TMM, RLE, or UPQ normalization method from edgeR package
#' @details Data is normalized and exported to wd for pandaPy()/lionessPy()
#' @return Returns clean, normalized expression data frame
#' @author Walid Nashashibi (\url{https://github.com/wnashash/})
#' @examples
#' \dontrun{
#' library(processPandaR)
#' normalized <- transform_norm(expression,"TMM")
#' }
#' @import edgeR
#' @import data.table
#' @import utils
#' @export
#'
transform_norm <- function(expData, method) {

  dgList <- edgeR::DGEList(counts = expData)

  # Normalize by TMM, RLE, or UPQ
  if(method == "TMM") {
    first <- edgeR::calcNormFactors(dgList, method = "TMM")
  } else if(method == "RLE") {
    first <- edgeR::calcNormFactors(dgList, method = "RLE")
  } else if(method == "UPQ") {
    first <- edgeR::calcNormFactors(dgList, method = "upperquartile")
  } else {
    stop("Method not supported. Please, choose between TMM, RLE, and UPQ.", call. = FALSE)
  }

  # Transform to normal distribution
  second <- log2(edgeR::cpm(first) + 1)
  second <- data.table::setDF(as.data.frame(second))

  # Write out results for pandaPy()/lionessPy()
  utils::write.table(second,
                     file=paste0(method, "_normalized.txt"),
                     row.names=TRUE, col.names=FALSE, sep="\t", quote=F)

  return(second)

}

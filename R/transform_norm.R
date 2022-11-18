#' @title Normalize Expression File
#' @description A simple function to normalize expression data and export
#' @param expData (Object) Name of expression dataframe object generated using readClean()
#' @param method (Character) Either TMM, RLE, or UPQ normalization method from edgeR package
#' @details Data is normalized and exported to working directory for pandaPy()/lionessPy()
#' @return Returns filtered, normalized expression dataframe
#' @author Walid Nashashibi (\url{https://github.com/wnashash/})
#' @examples
#' \dontrun{
#' library(processPandaR)
#' normalized <- transformNorm(expression,"TMM")
#' }
#' @import edgeR
#' @importFrom data.table setDF
#' @importFrom utils write.table
#' @export
#'
transform_norm <- function(expData, method) {

  dgList <- DGEList(counts = expData)

  # Normalize by TMM, RLE, or UPQ
  if(method == "TMM") {
    first <- calcNormFactors(dgList, method = "TMM")
  } else if(method == "RLE") {
    first <- calcNormFactors(dgList, method = "RLE")
  } else if(method == "UPQ") {
    first <- calcNormFactors(dgList, method = "upperquartile")
  } else {
    stop("Method not supported. Please, choose between TMM, RLE, and UPQ.", call. = FALSE)
  }

  # Transform to normal distribution
  second <- log2(cpm(first) + 1)
  second <- setDF(as.data.frame(second))

  # Write out results for pandaPy()/lionessPy()
  write.table(second,
              file=paste0(method, "_normalized.txt"),
              row.names=TRUE, col.names=FALSE, sep="\t", quote=F)

  return(second)
}

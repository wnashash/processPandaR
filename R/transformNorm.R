# Simple function to normalize expression matrix
transformNorm <- function(method) {

  ## First call readClean() function
  expDF <- readClean("AGM_expression_data.csv","csv","gene_name")

  ## Construct DGE list
  dgList <- DGEList(counts=expDF)

  ## Then normalize by TMM, RLE, or UPQ
  if(method == "TMM") {
    first <- calcNormFactors(dgList, method = "TMM")
  } else if(method == "RLE") {
    first <- calcNormFactors(dgList, method = "RLE")
  } else if(method == "UPQ") {
    first <- calcNormFactors(dgList, method = "upperquartile")
  }

  ## Transform to normal distribution
  second <- log2(cpm(first) + 1)

  return(second)
}

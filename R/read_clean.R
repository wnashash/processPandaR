#' @title Filter Expression File
#' @description A simple function to read and filter expression data
#' @param file (Character) File path, including file extension
#' @param geneName (Character) Name used for gene feature in expression data
#' @details Data is then normalized and exported to working directory, using transformNorm()
#' @return Returns filtered, raw expression dataframe
#' @author Walid Nashashibi (\url{https://github.com/wnashash/})
#' @examples
#' \dontrun{
#' library(processPandaR)
#' expression <- readClean("expression_data.csv","gene_name")
#' }
#' @import tools
#' @import genefilter
#' @importFrom data.table fread setDF
#' @export
#'
read_clean <- function(file, geneName){

  ext <- file_ext(file)

  if(ext == "rda") {
    expF <- load(file)
  } else {
    expF <- fread(file,
                  sep="auto",
                  header = T,
                  stringsAsFactors = F)
  }

  expData <- setDF(expF)

  # Filter rows with empty/duplicate gene name
  expData <- expData[!(is.na(expData[,geneName]) | expData[,geneName]==""), ]
  expData <- expData[!duplicated(expData[,geneName]), ]
  row.names(expData) <- expData[,geneName]
  expData <- expData[ ,!names(expData) %in% geneName]

  # Filter rows with low counts/variance
  thr <- 0.3
  expData <- expData[rowSums(expData) >= 10,]
  expData <- expData[rowSds(as.matrix(expData)) > thr*rowMeans(as.matrix(expData)),]

  return(expData)
}

#' @title Clean Gene Expression File
#' @description A simple function to read and clean gene expression data
#' @param file (Character) Gene expression file path, including file extension
#' @param geneName (Character) Name used for gene feature column in gene expression data
#' @return Returns clean gene expression data frame to global env and writes to wd
#' @author Walid Nashashibi (\url{https://github.com/wnashash/})
#' @examples
#' \dontrun{
#' library(processPandaR)
#' expression <- read_clean("expression.csv","gene")
#' }
#' @import tools
#' @import genefilter
#' @import data.table
#' @import utils
#' @export
#'
read_clean <- function(file, geneName){

  ext <- tools::file_ext(file)

  if(ext == "rda" | ext == "RData") {

    expF <- load(file)

  } else {

    expF <- data.table::fread(file,
                              sep="auto",
                              header=TRUE,
                              stringsAsFactors=FALSE)

  }

  expData <- data.table::setDF(expF)

  # Filter rows with empty/duplicate gene name
  expData <- expData[!(is.na(expData[,geneName]) | expData[,geneName]==""), ]
  expData <- expData[!duplicated(expData[,geneName]), ]
  row.names(expData) <- expData[,geneName]
  expData <- expData[ ,!names(expData) %in% geneName]

  # Filter rows with low counts/variance
  #thr <- 0.3
  #expData <- expData[rowSums(expData) >= 10,]
  #expData <- expData[genefilter::rowSds(as.matrix(expData)) > thr*rowMeans(as.matrix(expData)),]

  # Write out results for pandaPy()/lionessPy()
  utils::write.table(expData, file="expression_test.txt",
                     row.names=TRUE, col.names=FALSE, sep="\t", quote=F)

  return(expData)

}

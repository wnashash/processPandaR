#' @title Read Gene Expression File
#' @description A simple function to read and clean gene expression data
#' @param path (Character) Gene expression file path, including file extension
#' @param geneName (Character) Name used for gene feature column in gene expression data
#' @return Returns clean gene expression data frame to global env and writes to working directory
#' @author Walid Nashashibi (\url{https://github.com/wnashash/})
#' @examples
#' \dontrun{
#' library(processNetZoo)
#' exp_path <- system.file("extdata", "expression_test.csv", package = "processNetZoo", mustWork = TRUE)
#' expression <- read_gene(exp_path,'gene')
#' }
#' @import tools
#' @import data.table
#' @import utils
#' @export
#'
read_gene <- function(path, geneName){

  ext <- tools::file_ext(path)

  if(ext == "rda" | ext == "RData") {

    expF <- load(path)

  } else {

    expF <- data.table::fread(path,
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

  # Write out results for further processing
  utils::write.table(expData, file="expression.txt",
                     row.names=TRUE, col.names=FALSE, sep="\t", quote=F)

  return(expData)

}

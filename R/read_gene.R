#' @title Read Gene Expression File
#' @description A simple function to read and clean gene expression data
#' @param path (Character) Gene expression file path, including file extension
#' @param gene_column (Character) Name used for gene feature column in gene expression data
#' @param expr_file_out (Character) Optional argument to save output in working directory.  No
#' file will be saved when output_file is NULL.
#' @return Returns clean gene expression data frame to global env and writes to working directory
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
#' }
#' @import tools
#' @import data.table
#' @import utils
#' @export
#'
read_gene <- function(path, gene_column, expr_file_out = 'expression.txt'){

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
  expData <- expData[!(is.na(expData[,gene_column]) | expData[,gene_column]==""), ]
  expData <- expData[!duplicated(expData[,gene_column]), ]
  row.names(expData) <- expData[,gene_column]
  expData <- expData[ ,!names(expData) %in% gene_column]

  # Write out results for further processing
  if(!is.null(expr_file_out)){
    utils::write.table(expData, file = expr_file_out,
                     row.names=TRUE, col.names=FALSE, sep="\t", quote=F)
  }

  return(expData)

}

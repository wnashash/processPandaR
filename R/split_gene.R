#' @title Split Gene Expression File Using Metadata
#' @description A simple function to split gene expression data, using metadata treatment column
#' @param exp (Object) Name of expression data frame object generated using read_gene() or filter_gene()
#' @param metadata_file (Character) Metadata file path
#' @param id_column (Character) Metadata column with samples/IDs found in gene expression data
#' @param treatment_column (Character) Metadata column that designates control and treatment samples/IDs
#' @param expr_file_out (Character) Filename of the base group to write out results to.  No file will
#' be written when NULL.
#' @param comp_file_out (Character) Filename of the comparison group to write out results to.  No file
#' will be written when NULL.
#' @return Returns design and writes control (base) and treatment (comp) to working directory
#' @author Walid Nashashibi (\url{https://github.com/wnashash/})
#' @examples
#' \dontrun{
#' library(processNetZoo)
#'
#' exp_path <- system.file("extdata",
#'      "expression_test.csv",
#'      package = "processNetZoo",
#'      mustWork = TRUE
#'      )
#'
#' expression <- read_gene(exp_path,'gene')
#'
#' generate_histogram(expression,'both')
#'
#' filtered <- filter_gene(expression,'sum')
#'
#' meta_path <- system.file("extdata",
#'      "metadata.csv",
#'      package = "processNetZoo",
#'      mustWork = TRUE
#'      )
#'
#' design <- split_gene(filtered,meta_path,'SUBJECT_ID','RECURRENCE_ANY')
#' }
#' @import data.table
#' @import utils
#' @export
#'
split_gene <- function(exp,
                       metadata_file,
                       id_column,
                       treatment_column,
                       expr_file_out = 'expression_base.txt',
                       comp_file_out = 'expression_comp.txt') {

  metaF <- data.table::fread(metadata_file,
                             sep="auto",
                             header=TRUE,
                             stringsAsFactors=FALSE)

  metaData <- data.table::setDF(metaF)

  compFactors <- levels(as.factor(metaData[,treatment_column]))

  # Labels that identify baseline group samples/IDs
  baselineGrp <- which(metaData[,treatment_column] == compFactors[1])
  baseGrp     <- metaData[baselineGrp,]

  geneExpBase <- exp[,baseGrp[,id_column]]

  if(!is.null(expr_file_out)){
    utils::write.table(geneExpBase, file = expr_file_out,
                     row.names=TRUE, col.names=FALSE, sep="\t", quote=FALSE)
  }

  # Labels that identify the comparison group samples/IDs
  compareGrp  <- which(metaData[,treatment_column] == compFactors[2])
  compGrp     <- metaData[compareGrp,]

  geneExpComp <- exp[,compGrp[,id_column]]

  if(!is.null(comp_file_out)){
    utils::write.table(geneExpComp, file = comp_file_out,
                     row.names=TRUE, col.names=FALSE, sep="\t", quote=FALSE)
  }

  # Generate design for monster algorithm
  metaData <- metaData[order(match(metaData[,id_column],colnames(exp))),]
  design <- ifelse(metaData[,treatment_column] == compFactors[1], 0, 1)

  return(design)

}


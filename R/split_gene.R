#' @title Split Gene Expression File Using Metadata
#' @description A simple function to split gene expression data, using metadata treatment column
#' @param geneExp (Object) Name of expression data frame object generated using read_gene() or filter_gene()
#' @param metaPath (Character) Metadata file path
#' @param simpleID (Character) Metadata column with samples/IDs found in gene expression data
#' @param treatment (Character) Metadata column that designates control and treatment samples/IDs
#' @return Returns design and writes control (base) and treatment (comp) to working directory
#' @author Walid Nashashibi (\url{https://github.com/wnashash/})
#' @examples
#' \dontrun{
#' library(processNetZoo)
#' expression <- read_gene('extdata/expression_test.csv','gene')
#' generate_histogram(expression,'both')
#' filtered <- filter_gene(expression,'sum')
#' design <- split_gene(filtered,'extdata/metadata.csv','SUBJECT_ID','RECURRENCE_ANY')
#' }
#' @import data.table
#' @import utils
#' @export
#'
split_gene <- function(geneExp,metaPath,simpleID,treatment) {

  metaF <- data.table::fread(metaPath,
                             sep="auto",
                             header=TRUE,
                             stringsAsFactors=FALSE)

  metaData <- data.table::setDF(metaF)

  compFactors <- levels(as.factor(metaData[,treatment]))

  # Labels that identify baseline group samples/IDs
  baselineGrp <- which(metaData[,treatment] == compFactors[1])
  baseGrp     <- metaData[baselineGrp,]

  geneExpBase <- geneExp[,baseGrp[,simpleID]]

  utils::write.table(geneExpBase, file="expression_base.txt",
                     row.names=TRUE, col.names=FALSE, sep="\t", quote=FALSE)

  # Labels that identify the comparison group samples/IDs
  compareGrp  <- which(metaData[,treatment] == compFactors[2])
  compGrp     <- metaData[compareGrp,]

  geneExpComp <- geneExp[,compGrp[,simpleID]]

  utils::write.table(geneExpComp, file="expression_comp.txt",
                     row.names=TRUE, col.names=FALSE, sep="\t", quote=FALSE)

  # Generate design for monster algorithm
  metaData <- metaData[order(match(metaData[,simpleID],colnames(geneExp))),]
  design <- ifelse(metaData[,treatment] == compFactors[1], 0, 1)

  return(design)

}


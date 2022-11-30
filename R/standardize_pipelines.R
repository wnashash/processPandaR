#' @title File Prep and Run Pipeline
#' @description A simple function to prepare data and run pipeline of choice
#' @param pipeline (Character) Choose between panda, lioness, pandaCondor, and pandaAlpaca
#' @param fileGeneExp (Character) File path of gene expression file
#' @param colGeneExp (Character) Name of gene feature column in expression file
#' @param fileMetadata (Character) File path of metadata file
#' @param colMetadata (Character) Name of treatment feature column in metadata file
#' @return Returns output of selected pipeline to be used for visualization
#' @author Walid Nashashibi (\url{https://github.com/wnashash/})
#' @examples
#' \dontrun{
#' library(processPandaR)
#' alpacaResult <- run_pipeline("pandaAlpaca","expression_test.csv","gene","metadata.csv","RECURRENCE_ANY")
#' }
#' @import tools
#' @import data.table
#' @import utils
#' @export
#'
run_pipeline <- function(pipeline,fileGeneExp,colGeneExp,fileMetadata=NULL,colMetadata=NULL){

  ext <- tools::file_ext(fileGeneExp)

  if(ext == "rda" | ext == "RData") {

    expF <- load(fileGeneExp)

  } else {

    expF <- data.table::fread(fileGeneExp,
                              sep="auto",
                              header=TRUE,
                              stringsAsFactors=FALSE)

  }

  expData <- data.table::setDF(expF)

  # Filter rows with empty/duplicate gene name
  expData <- expData[!(is.na(expData[,colGeneExp]) | expData[,colGeneExp]==""), ]
  expData <- expData[!duplicated(expData[,colGeneExp]), ]
  row.names(expData) <- expData[,colGeneExp]
  expData <- expData[ ,!names(expData) %in% colGeneExp]

  # Filter rows with low counts/variance
  #thr <- 0.3
  #expData <- expData[rowSums(expData) >= 10,]
  #expData <- expData[genefilter::rowSds(as.matrix(expData)) > thr*rowMeans(as.matrix(expData)),]

  # Write out results for pandaPy()/lionessPy()
  utils::write.table(expData, file="expression_test.txt",
                     row.names=TRUE, col.names=FALSE, sep="\t", quote=FALSE)

  # Run pandaPy() for use in pipelines - not wrapping in with statement for now
  resultPanda <- netZooR::pandaPy("expression_test.txt",
                                  "OV_Motif.txt",
                                  "OV_PPI.txt",
                                  save_tmp = FALSE,
                                  modeProcess = "intersection")
  if(pipeline == "panda") {

    result <- resultPanda$panda

  } else if(pipeline == "lioness") {

    withr::with_file("expression_test.txt",
                     {resultLion <- netZooR::lionessPy("expression_test.txt",
                                                       "OV_Motif.txt",
                                                       "OV_PPI.txt",
                                                       save_tmp = FALSE,
                                                       modeProcess = "intersection")})

    result <- resultLion[c(3:ncol(resultLion))]
    colnames(result) <- colnames(expData)
    result <- cbind(resultLion[,c("tf","gene")], result)

  } else if(pipeline == "pandaCondor") {

    resultCondor <- netZooR::pandaToCondorObject(resultPanda$panda)

    # Not sure if this is the result that should be returned for viz
    result <- netZooR::condorCluster(resultCondor)

  } else if(pipeline == "pandaAlpaca") {

    metaF <- data.table::fread(fileMetadata,
                               sep="auto",
                               header=TRUE,
                               stringsAsFactors=FALSE)

    metaData <- data.table::setDF(metaF)

    compFactors <- levels(as.factor(metaData[,colMetadata]))

    # Labels (row numbers) that identify baseline group patients
    baselineGrp <- which(metaData[,colMetadata] == compFactors[1])
    baseGrp     <- metaData[baselineGrp,]
    geneExpBase <- expData[,baseGrp$SUBJECT_ID] # need to standardize this

    utils::write.table(geneExpBase, file="expression_base.txt",
                       row.names=TRUE, col.names=FALSE, sep="\t", quote=FALSE)

    withr::with_file("expression_base.txt",
                     {controlPanda <- netZooR::pandaPy("expression_base.txt",
                                                       "OV_Motif.txt",
                                                       "OV_PPI.txt",
                                                       save_tmp = FALSE,
                                                       modeProcess = "intersection")$panda})

    # Labels (row numbers) that identify the comparison group patients
    compareGrp  <- which(metaData[,colMetadata] == compFactors[2])
    compGrp     <- metaData[compareGrp,]
    geneExpComp <- expData[,compGrp$SUBJECT_ID] # need to standardize this

    utils::write.table(geneExpComp, file="expression_comp.txt",
                       row.names=TRUE, col.names=FALSE, sep="\t", quote=FALSE)

    withr::with_file("expression_comp.txt",
                     {treatPanda <- netZooR::pandaPy("expression_comp.txt",
                                                     "OV_Motif.txt",
                                                     "OV_PPI.txt",
                                                     save_tmp = FALSE,
                                                     modeProcess = "intersection")$panda})

    result <- netZooR::pandaToAlpaca(treatPanda, controlPanda)

  } else {

    stop("Pipeline not supported. Please, choose between panda, lioness, pandaCondor, and pandaAlpaca.", call. = FALSE)

  }

  return(result)

}

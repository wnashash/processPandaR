#' @title NetZooR Pipeline Processing Steps
#' @description A series of simple functions to prepare data and run algorithm of choice
#' @param pipeline (Character) Choose between panda, lioness, condor, alpaca, and monster
#' @param pathGene (Character) File path of gene expression file
#' @param geneName (Character) Name of gene feature column in expression file
#' @param pathMeta (Character) File path of metadata file for alpaca and monster
#' @param simpleID (Character) Name of sample ID feature column in metadata file
#' @param treatment (Character) Name of treatment feature column in metadata file
#' @return Returns output of selected algorithm to be used for visualization or second algorithm
#' @author Walid Nashashibi (\url{https://github.com/wnashash/})
#' @examples
#' \dontrun{
#' library(processPandaR)
#' resultPanda <- algo_select("panda",exp,"expression_test.csv")
#' }
#' @import tools
#' @import data.table
#' @import utils
#' @export
#'
read_geneExp <- function(pathGene, geneName){

  ext <- tools::file_ext(pathGene)

  if(ext == "rda" | ext == "RData") {

    expF <- load(pathGene)

  } else {

    expF <- data.table::fread(pathGene,
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

  # Write out results for pandaPy()/lionessPy()
  utils::write.table(expData, file="expression_test.txt",
                     row.names=TRUE, col.names=FALSE, sep="\t", quote=F)

  return(expData)

}



hist_geneExp <- function(expData,typeFilt) {

  if(typeFilt=='sum') {
    hist(rowSums(expData))
  } else if(typeFilt=='var') {
    hist(apply(expData,1,var))
  } else {
    stop("Filter not supported. Please, choose between sum and var.", call. = FALSE)
  }

}



filt_geneExp <- function(expData,typeFilt,thresh=NULL) {

  if(typeFilt=='sum') {

    expData <- expData[rowSums(expData) >= 10,]

  } else if(typeFilt=='var') {

    if(is.null(thresh)) {
      thr <- 0.3
    } else {
      thr <- thresh
    }

    expData <- expData[genefilter::rowSds(as.matrix(expData)) > thr*rowMeans(as.matrix(expData)),]

  } else if(typeFilt=='both') {

    expData <- expData[rowSums(expData) >= 10,]

    if(is.null(thresh)) {
      thr <- 0.3
    } else {
      thr <- thresh
    }

    expData <- expData[genefilter::rowSds(as.matrix(expData)) > thr*rowMeans(as.matrix(expData)),]

  } else {

    stop("Filter not supported. Please, choose between sum, var, or both.", call. = FALSE)

  }

  # Write out results for pandaPy()/lionessPy()
  utils::write.table(expData, file="expression_test.txt",
                     row.names=TRUE, col.names=FALSE, sep="\t", quote=F)

  return(expData)

}



use_metaData <- function(expData,pathMeta,simpleID,treatment) {

  metaF <- data.table::fread(pathMeta,
                             sep="auto",
                             header=TRUE,
                             stringsAsFactors=FALSE)

  metaData <- data.table::setDF(metaF)

  compFactors <- levels(as.factor(metaData[,treatment]))

  # Labels (row numbers) that identify baseline group patients
  baselineGrp <- which(metaData[,treatment] == compFactors[1])
  baseGrp     <- metaData[baselineGrp,]
  geneExpBase <- expData[,baseGrp[,simpleID]]

  utils::write.table(geneExpBase, file="expression_base.txt",
                     row.names=TRUE, col.names=FALSE, sep="\t", quote=FALSE)

  # Labels (row numbers) that identify the comparison group patients
  compareGrp  <- which(metaData[,treatment] == compFactors[2])
  compGrp     <- metaData[compareGrp,]
  geneExpComp <- expData[,compGrp[,simpleID]] # need to standardize this

  utils::write.table(geneExpComp, file="expression_comp.txt",
                     row.names=TRUE, col.names=FALSE, sep="\t", quote=FALSE)

}



algo_select <- function(algorithm,expData,basePath,compPath=NULL) {

  if(algorithm=='panda') {

    withr::with_file("expression_test.txt",
                     {resultPanda <- netZooR::pandaPy("expression_test.txt",
                                                      "OV_Motif.txt",
                                                      "OV_PPI.txt",
                                                      save_tmp = FALSE,
                                                      modeProcess = "intersection")$panda})

    result <- resultPanda

  } else if(algorithm=='lioness') {

    withr::with_file("expression_test.txt",
                     {resultLion <- netZooR::lionessPy("expression_test.txt",
                                                       "OV_Motif.txt",
                                                       "OV_PPI.txt",
                                                       save_tmp = FALSE,
                                                       modeProcess = "intersection")})

    result <- resultLion[c(3:ncol(resultLion))]
    colnames(result) <- colnames(expData)
    result <- cbind(resultLion[,c("tf","gene")], result)

  } else if(algorithm=='condor') {

    withr::with_file("expression_test.txt",
                     {resultPanda <- netZooR::pandaPy("expression_test.txt",
                                                      "OV_Motif.txt",
                                                      "OV_PPI.txt",
                                                      save_tmp = FALSE,
                                                      modeProcess = "intersection")$panda})

    resultCondor <- netZooR::pandaToCondorObject(resultPanda)
    resultCondor <- netZooR::condorCluster(resultCondor)

    result <- resultCondor

  } else if(algorithm=='alpaca') {

    withr::with_file("expression_base.txt",
                     {controlPanda <- netZooR::pandaPy("expression_base.txt",
                                                       "OV_Motif.txt",
                                                       "OV_PPI.txt",
                                                       save_tmp = FALSE,
                                                       modeProcess = "intersection")$panda})

    withr::with_file("expression_comp.txt",
                     {treatPanda <- netZooR::pandaPy("expression_comp.txt",
                                                     "OV_Motif.txt",
                                                     "OV_PPI.txt",
                                                     save_tmp = FALSE,
                                                     modeProcess = "intersection")$panda})

    tableNet <- merge(controlPanda[,-3],treatPanda[,-3],by = c("TF","Gene"))

    resultAlpaca <- netZooR::alpaca(tableNet,'./alpaca')
    resultCrane  <- netZooR::alpacaCrane(tableNet,resultAlpaca)

    resultTF <- resultCrane$TF
    resultG  <- resultCrane$Gene

    result <- list() # placeholder

  } else if(algorithm=='monster') {

    result <- list() # placeholder

  }

  return(result)

}

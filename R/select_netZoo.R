#' @title Select NetZooR Algorithm
#' @description A simple function to run NetZooR algorithm of choice
#' @param algorithm (Character) Choose between panda, lioness, condor, alpaca, and monster
#' @param geneExp (Object) Name of expression data frame object generated using read_gene() or filter_gene()
#' @param basePath (Character) Gene expression base group file path
#' @param compPath (Character) Gene expression comparison group file path
#' @param startSample (Numeric) Value to indicate lionessPy first sample - default is 1
#' @param endSample (Numeric) Value to indicate lionessPy last sample - default is 'None'
#' @details In order to run alpaca or monster, split_gene() must be run first
#' @return Returns output of selected algorithm to be used for visualization or second algorithm
#' @author Walid Nashashibi (\url{https://github.com/wnashash/})
#' @examples
#' \dontrun{
#' library(processNetZoo)
#' exp_path <- system.file("extdata", "expression_test.csv", package = "processNetZoo", mustWork = TRUE)
#' expression <- read_gene(exp_path,'gene')
#' result <- select_NetZoo('lionessPy',expression,'expression.txt',startSample=11,endSample=20)
#' }
#' @import netZooR
#' @import data.table
#' @export
#'
select_NetZoo <- function(algorithm,geneExp,basePath,compPath=NULL,startSample=1,endSample='None') {

  motif_path <- system.file("extdata", "OV_Motif.txt", package = "processNetZoo", mustWork = TRUE)
  ppi_path <- system.file("extdata", "OV_PPI.txt", package = "processNetZoo", mustWork = TRUE)

  if(algorithm == 'pandaPy') {

    resultPanda <- netZooR::pandaPy(basePath,
                                    motif_path,
                                    ppi_path,
                                    save_tmp = FALSE,
                                    modeProcess = "intersection")$panda

    result <- resultPanda

  } else if(algorithm == 'lionessPy') {

    resultLion <- netZooR::lionessPy(basePath,
                                     motif_path,
                                     ppi_path,
                                     save_tmp = FALSE,
                                     modeProcess = "intersection",
                                     start_sample = startSample,
                                     end_sample = endSample)

    result <- resultLion[c(3:ncol(resultLion))]

    if(endSample == 'None') {
      geneTruncated <- geneExp[c(startSample:ncol(geneExp))]
    } else {
      geneTruncated <- geneExp[c(startSample:endSample)]
    }

    colnames(result) <- colnames(geneTruncated)
    result <- cbind(resultLion[,c("tf","gene")], result)

    data.table::setnames(result, old = c('tf','gene'), new = c('TF','Gene'))

  } else if(algorithm == 'condor') {

    resultPanda <- netZooR::pandaPy(basePath,
                                    motif_path,
                                    ppi_path,
                                    save_tmp = FALSE,
                                    modeProcess = "intersection")$panda

    resultCondor <- netZooR::pandaToCondorObject(resultPanda)
    resultCondor <- netZooR::condorCluster(resultCondor)
    result <- resultCondor

  } else if(algorithm == 'alpaca') {

    controlPanda <- netZooR::pandaPy(basePath,
                                     motif_path,
                                     ppi_path,
                                     save_tmp = FALSE,
                                     modeProcess = "intersection")$panda

    treatPanda   <- netZooR::pandaPy(compPath,
                                     motif_path,
                                     ppi_path,
                                     save_tmp = FALSE,
                                     modeProcess = "intersection")$panda

    tableNet <- merge(controlPanda[,-3],treatPanda[,-3],by = c("TF","Gene"))

    resultAlpaca <- netZooR::alpaca(tableNet,'./alpaca')
    resultCrane  <- netZooR::alpacaCrane(tableNet,resultAlpaca)
    result <- resultCrane

  } else if(algorithm == 'monster') {

    motif <- data.table::fread(motif_path,
                               sep="auto",
                               header=FALSE,
                               stringsAsFactors=FALSE)

    motif <- data.table::setDF(motif)

    monsterResult <- netZooR::monster(geneExp,
                                      design,
                                      motif,
                                      nullPerms = 10,
                                      mode = 'buildNet') # different set up for regNet

    result <- monsterResult

  }

  return(result)

}

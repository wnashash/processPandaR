#' @title Select NetZooR Algorithm
#' @description A simple function to run NetZooR algorithm of choice
#' @param alg (Character) Choose between panda, lioness, condor, alpaca, and monster
#' @param exp (Object) Name of expression data frame object generated using read_gene() or filt_gene()
#' @param base (Character) Gene expression base group file path
#' @param comp (Character) Gene expression comparison group file path
#' @details In order to run alpaca or monster, split_gene() must be run first
#' @return Returns output of selected algorithm to be used for visualization or second algorithm
#' @author Walid Nashashibi (\url{https://github.com/wnashash/})
#' @examples
#' \dontrun{
#' library(processNetZoo)
#' expression <- read_gene('expression_test.csv','gene')
#' design <- split_gene(expression,'metadata.csv','SUBJECT_ID','RECURRENCE_ANY')
#' result <- select_NetZoo('alpaca',expression,'expression_base.txt','expression_comp.txt')
#' }
#' @import netZooR
#' @import data.table
#' @export
#'
select_NetZoo <- function(alg,exp,base,comp=NULL) {

  if(alg=='panda') {

    resultPanda <- netZooR::pandaPy(base,
                                    "OV_Motif.txt",
                                    "OV_PPI.txt",
                                    save_tmp = FALSE,
                                    modeProcess = "intersection")$panda

    result <- resultPanda

  } else if(alg=='lioness') {

    resultLion <- netZooR::lionessPy(base,
                                     "OV_Motif.txt",
                                     "OV_PPI.txt",
                                     save_tmp = FALSE,
                                     modeProcess = "intersection")

    result <- resultLion[c(3:ncol(resultLion))]
    colnames(result) <- colnames(exp)
    result <- cbind(resultLion[,c("tf","gene")], result)

  } else if(alg=='condor') {

    resultPanda <- netZooR::pandaPy(base,
                                    "OV_Motif.txt",
                                    "OV_PPI.txt",
                                    save_tmp = FALSE,
                                    modeProcess = "intersection")$panda

    resultCondor <- netZooR::pandaToCondorObject(resultPanda)
    resultCondor <- netZooR::condorCluster(resultCondor)
    result <- resultCondor

  } else if(alg=='alpaca') {

    controlPanda <- netZooR::pandaPy(base,
                                     "OV_Motif.txt",
                                     "OV_PPI.txt",
                                     save_tmp = FALSE,
                                     modeProcess = "intersection")$panda

    treatPanda   <- netZooR::pandaPy(comp,
                                     "OV_Motif.txt",
                                     "OV_PPI.txt",
                                     save_tmp = FALSE,
                                     modeProcess = "intersection")$panda

    tableNet <- merge(controlPanda[,-3],treatPanda[,-3],by = c("TF","Gene"))

    resultAlpaca <- netZooR::alpaca(tableNet,'./alpaca')
    resultCrane  <- netZooR::alpacaCrane(tableNet,resultAlpaca)
    result <- resultCrane

  } else if(alg=='monster') {

    motif <- data.table::fread("OV_Motif.txt",
                               sep="auto",
                               header=FALSE,
                               stringsAsFactors=FALSE)

    motif <- data.table::setDF(motif)

    monsterResult <- netZooR::monster(exp,
                                      design,
                                      motif,
                                      nullPerms = 10,
                                      mode = 'buildNet') # different set up for regNet

    result <- monsterResult

  }

  return(result)

}


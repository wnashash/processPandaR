#' @title Select NetZooR Algorithm
#' @description A simple function to run NetZooR algorithm of choice
#' @param algorithm (Character) Choose between 'pandaPy', 'lionessPy', 'condor', 'alpaca', and 'monster'
#' @param expr_file (Character) Gene expression base group filepath
#' @param expr (Object) Name of expression data.frame object generated using `read_gene()` or `filter_gene()`
#' @param motif_file (Character) value indicating which motif data to use.  Options include NULL for Pearson correlation
#' matrix, 'default' for built-in data, or a user-specified filepath.
#' @param ppi_file (Character) value indicating which protein-protein interaction data to use.  Options include
#' NULL for no additional interaction, 'default' for built-in data, or a user-specified filepath.
#' @param comp_file (Character) Gene expression comparison group filepath for 'alpaca' networks
#' @param start_sample (Numeric) Value to indicate 'lionessPy' first sample - default is 1
#' @param end_sample (Numeric) Value to indicate 'lionessPy' last sample - default is 'None'
#' @param design (Numeric) Binary vector for 'monster' algorithm case control partition.
#' @details In order to run alpaca or monster, split_gene() must be run first
#' @return Returns output of selected algorithm to be used for visualization or second algorithm
#' @author Walid Nashashibi (\url{https://github.com/wnashash/})
#' @examples
#' \dontrun{
#' library(processNetZoo)
#'
#' # -- pandaPy --
#' exp_path <- system.file("extdata",
#'      "expr4_matched.txt",
#'      package = "netZooR",
#'      mustWork = TRUE)
#'
#' result <- select_NetZoo(
#'      algorithm   = 'pandaPy',
#'      expr_file   = exp_path,
#'      motif_file  = NULL,
#'      ppi_file    = NULL
#' )
#'
#'
#' # -- lionessPy --
#' exp_path <- system.file("extdata",
#'      "expression_test.csv",
#'      package = "processNetZoo",
#'      mustWork = TRUE)
#'
#' expression <- read_gene(exp_path,'gene')
#'
#' result <- select_NetZoo(
#'      algorithm    = 'lionessPy',
#'      expr         = expression,
#'      expr_file    = 'expression.txt',
#'      start_sample = 11,
#'      end_sample   = 20
#')
#' }
#' @import netZooR
#' @import data.table
#' @export
#'
select_NetZoo <- function(algorithm,
                          expr_file    = NULL,
                          expr         = NULL,
                          motif_file   = 'default',
                          ppi_file     = 'default',
                          comp_file    = NULL,
                          start_sample = 1,
                          end_sample   = 'None',
                          design       = NULL) {

  # set up motif and PPI files
  if(!is.null(motif_file) && motif_file == 'default'){
    motif_file <- system.file("extdata", "OV_Motif.txt", package = "processNetZoo", mustWork = TRUE)
  }

  if(!is.null(ppi_file) && ppi_file == 'default'){
    ppi_file <- system.file("extdata", "OV_PPI.txt", package = "processNetZoo", mustWork = TRUE)
  }

  if(algorithm == 'pandaPy') {

    resultPanda <- netZooR::pandaPy(expr_file   = expr_file,
                                    motif_file  = motif_file,
                                    ppi_file    = ppi_file,
                                    save_tmp    = FALSE,
                                    modeProcess = "intersection")$panda

    result <- resultPanda

  } else if(algorithm == 'lionessPy') {

    resultLion <- netZooR::lionessPy(expr_file    = expr_file,
                                     motif_file   = motif_file,
                                     ppi_file     = ppi_file,
                                     save_tmp     = FALSE,
                                     modeProcess  = "intersection",
                                     start_sample = start_sample,
                                     end_sample   = end_sample,
                                     save_single_network = FALSE,
                                     save_dir = NULL)

    result <- resultLion[c(3:ncol(resultLion))]

    if(end_sample == 'None') {
      geneTruncated <- expr[c(start_sample:ncol(expr))]
    } else {
      geneTruncated <- expr[c(start_sample:end_sample)]
    }

    colnames(result) <- colnames(geneTruncated)
    result <- cbind(resultLion[,c("tf","gene")], result)

    data.table::setnames(result, old = c('tf','gene'), new = c('TF','Gene'))

  } else if(algorithm == 'condor') {

    resultPanda <- netZooR::pandaPy(expr_file   = expr_file,
                                    motif_file  = motif_file,
                                    ppi_file    = ppi_file,
                                    save_tmp    = FALSE,
                                    modeProcess = "intersection")$panda

    resultCondor <- netZooR::pandaToCondorObject(resultPanda)
    resultCondor <- netZooR::condorCluster(resultCondor)
    result <- resultCondor

  } else if(algorithm == 'alpaca') {

    controlPanda <- netZooR::pandaPy(expr_file   = expr_file,
                                     motif_file  = motif_file,
                                     ppi_file    = ppi_file,
                                     save_tmp    = FALSE,
                                     modeProcess = "intersection")$panda

    treatPanda   <- netZooR::pandaPy(expr_file   = comp_file,
                                     motif_file  = motif_file,
                                     ppi_file    = ppi_file,
                                     save_tmp    = FALSE,
                                     modeProcess = "intersection")$panda

    tableNet <- merge(controlPanda[,-3],treatPanda[,-3],by = c("TF","Gene"))

    resultAlpaca <- netZooR::alpaca(tableNet,'./alpaca')
    resultCrane  <- netZooR::alpacaCrane(tableNet,resultAlpaca)
    result <- resultCrane

  } else if(algorithm == 'monster') {

    motif <- data.table::fread(motif_file,
                               sep="auto",
                               header=FALSE,
                               stringsAsFactors=FALSE)

    motif <- data.table::setDF(motif)

    monsterResult <- netZooR::monster(expr,
                                      design,
                                      motif,
                                      nullPerms = 10,
                                      mode = 'buildNet') # different set up for regNet

    result <- monsterResult

  }

  return(result)

}

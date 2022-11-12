# Simple function to read and filter expression data
readClean <- function(file, format, geneName){

  ## Import file formats rda, csv, tsv, or txt
  if(format == "rda") {
    exp <- load(file = file)
  } else {
    exp <- fread(file = file,
                 sep="auto",
                 header = T,
                 stringsAsFactors = F)
  }

  ## Create DF - maybe unnecessary
  expDF <- setDF(exp)

  ## Filter rows with empty/duplicate gene name
  expDF <- expDF[!(is.na(expDF[,geneName]) | expDF[,geneName]==""), ]
  expDF <- expDF[!duplicated(expDF[,geneName]), ]
  row.names(expDF) <- expDF[,geneName]
  expDF <- expDF[ ,!names(expDF) %in% geneName]

  ## Filter rows with low counts/variance
  thr <- 0.3
  expDF <- expDF[rowSums(expDF) >= 10,]
  expDF <- expDF[rowSds(as.matrix(expDF)) > thr*rowMeans(as.matrix(expDF)),]

  return(expDF)
}

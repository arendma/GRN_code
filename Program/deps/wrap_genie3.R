wrap_genie3 <- function (dat, parallel, tfs=FALSE) {
  #takes a matrix of gene expression values and calculates genie 3 network for all gene combinations
  #gene names must be specified in rownames(data)
  library(GENIE3)
  if (tfs!=FALSE) {
    if (parallel>1) {
      weightMat <- GENIE3(dat, regulators = tfs,  nCores=parallel, verbose=TRUE)
    }
    else if (parallel==1) {
      weightMat <- GENIE3(dat,regulators =tfs, verbose=TRUE)
    }
    else {stop('parallel must be a positive integer')}
  }
  else {
    if (parallel>1) {
      weightMat <- GENIE3(dat, nCores=parallel, verbose=TRUE)
    }
    else if (parallel==1) {
      weightMat <- GENIE3(dat, verbose=TRUE)
    }
    else {stop('parallel must be a positive integer')}
  }
  linkList <- getLinkList(weightMat)
  genie_all <- linkList
  
  
  colnames(genie_all) <- c("from", "to", "weight")
  
  genie_all$from <- as.character(genie_all$from)
  genie_all$to <- as.character(genie_all$to)
  genie_all$weight <- as.numeric(genie_all$weight)
  
  return(genie_all)
}
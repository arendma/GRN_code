wrap_gennet <- function(dat, tfs=FALSE) {
  #takes a matrix of gene expression values and calculates network for all gene combinations based on ggm estimates
  #gene names must be specified in rownames(data) if a list of transcription factors (same ID as rownames) is supplied
  #the generated network is filtered for edges starting from these
  library(GeneNet)
  inferred.pcor <- ggm.estimate.pcor(t(dat))
  
  test.results <- network.test.edges(inferred.pcor,fdr=F,direct=F,plot=F)
  
  net <- test.results 
  # net <- extract.network(test.results,method.ggm="qval",cutoff.ggm=0.05)
  
  ggm_list <- data.frame(from=rownames(inferred.pcor)[net$node1],
                         to=rownames(inferred.pcor)[net$node2],
                         weight=net$pcor)
  if (length(tfs)>1 && class(tfs)=="character") {
    ggm_list <- ggm_list[ggm_list$from %in% tfs,]
  } else if(length(tfs)>1 || tfs!=FALSE) {stop("tfs argument is either FALSE or a character vector of dat row indexes")}
  
  #Taken from nooshins inferrence code but I don't now what g_rof_syn is (undeclared variable)
  #ggm_list$from <- tolower(as.character(g_orf_syn[ggm_list$from]))
  #ggm_list$to <- tolower(as.character(g_orf_syn[ggm_list$to]))
  #head(ggm_list)
  return(ggm_list)
}
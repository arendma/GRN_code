wrap_aracne <- function(dat, tfs=FALSE) {
  #takes a matrix of gene expression values and calculates network for all gene combinations based on aracne algotihm
  #and spearman estimates
  #gene names must be specified in rownames(data), if only interactions between TFs and Genes should be considered supply
  # a character vector of TF gene names as second argument
library(minet)# ARACNE and CLR
library(igraph)
mi_all <- build.mim(t(dat),estimator="spearman")
#mi_all <- mi_all[genes_to_use,genes_to_use]
net_all <- aracne(mi_all)
if (tfs!=FALSE) {
  net_all[!(rownames(net_all)%in%tfs),] <- 0
}

g_ara <- graph.adjacency(net_all,weighted=T)
a_edge <- get.data.frame(g_ara)
a_edge <- a_edge[order(a_edge$weight,decreasing=T),]
#Taken from nooshins inferrence code but I don't now what g_rof_syn is (undeclared variable)
#a_edge$to <- tolower(as.character(g_orf_syn[a_edge$to]))
#a_edge$from <- tolower(as.character(g_orf_syn[a_edge$from]))
return(a_edge)
}
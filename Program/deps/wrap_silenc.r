wrap_silenc <- function(dat, tfs=FALSE) {
  library(amap)
  source("deps/preprocess_silencing.r")
  source("deps/silencing.r")
  cor_ds <- cor(t(dat), method='pearson')
  diag(cor_ds) <- 1
  aa <- matrix(rnorm(length(cor_ds),mean=0,sd=1),ncol=ncol(cor_ds),nrow=nrow(cor_ds))/100
  
  if (tfs!=FALSE) {
    S <- Spreprocess(cor_ds+aa, tfs=NULL, nontfs=rownames(cor_ds)[!(rownames(cor_ds) %in% tfs)])
  } else {
    S <- Spreprocess(cor_ds+aa,tfs=NULL,nontfs=NULL) # variable tfs was unused in current code made nontfs also unused since I didn't understand the usage see preprocess_silencing.r
  }
  s<-S
  colnames(s) <- rownames(s)
  print(all(s[!(rownames(s) %in% tfs),]==0))
  if (tfs!=FALSE) { # unnescessary since this is already done in the preprocess silencing step
  s[!(rownames(s) %in% tfs),] <- 0
  }
  g_S <- graph.adjacency(s,weighted=T)
  s_edge <- get.data.frame(g_S)
  s_edge <- s_edge[order(s_edge$weight,decreasing=T),]
  #Taken from nooshins inferrence code but I don't now what g_rof_syn is (undeclared variable)
  #s_edge$from <- tolower(as.character(g_orf_syn[s_edge$from]))
  #s_edge$to <- tolower(as.character(g_orf_syn[s_edge$to]))
  return(s_edge)  
}
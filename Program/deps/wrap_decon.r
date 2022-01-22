wrap_decon <- function(dat,silent,tfs=FALSE) {
  #takes a expression data and c

library(amap)
source("deps/preprocess_deconvolution.r")
source("deps/deconvolution.r")
cor_ds <- cor(t(dat), method='pearson')
D <- Dpreprocess(cor_ds,NULL,NULL, silent) # second two arguments are not used in current version of nooshins code
d<-D
if (tfs!=FALSE) {
d[rownames(dat)[!(rownames(dat)%in% tfs)],] <- 0
}
g_D <- graph.adjacency(d,weighted=T)
d_edge <- get.data.frame(g_D)
#Taken from nooshins inferrence code but I don't now what g_rof_syn is (undeclared variable)
#from <- tolower(as.character(g_orf_syn[d_edge$from]))
#d_edge$to <- tolower(as.character(g_orf_syn[d_edge$to]))
return(d_edge)
}
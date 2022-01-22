## Script to create the elastic net derived grn
source('deps/wrap_elnetv3.r')
library('edgeR')
set.seed(2019)
options(stringsAsFactors = FALSE)
#load replicate wise counts tMM and voom standaradised merch1519arma17trd HTSeq have been filtered for genes with cpm>1 in at least 9 reps
load('../Data/20191112me1915arma17_exp.obj')


gene_dat <- t(scale(t(me1519arma17_exp)))
#check if transformation has introduced NA values
if(sum(is.na(gene_dat)!=0)) {stop('transformation introduced NA values check expression data')}

#Collapse to the median
getmedian <- function(x, gene_dat){ # ADJUSTED FOR REPLICATE INDEX STARTING WITH _X
  #check if data has replicas 
  if (sum(gsub('_[[:digit:]]', '', colnames(gene_dat))==x)>1) {
    return(apply(gene_dat[, gsub('_[[:digit:]]', '', colnames(gene_dat))==x],1, median))
  } else {return(gene_dat[, gsub('_[[:digit:]]', '', colnames(gene_dat))==x])}
} 
gene_dat_median <- sapply(unique(gsub('_[[:digit:]]', '', colnames(gene_dat))), getmedian, gene_dat=gene_dat)

#import TF annotation - names should be the same as rownames of the read data
TF <- read.delim('../Data/BLAST_curated_JIN_TF.txt', header =TRUE)
tfs <- TF$final_ID
tfs <- tfs[tfs %in% rownames(gene_dat)]

#Do elastic net regression analysis for each TF using all other TFs as predictors 
#THIS STEP TAKES VERY LONG - RUN ON SERVER OR wITH MORE THEN 8 THREADS
elnet_res <- wrap_elnet(gene_dat_median, resdir='../Results/Bignet/elnet/', thrsh=0, tfs=tfs,parallel=8)
elnet_all <- elnet_res[,c('Gene.ID', 'predicted', 'rel.coeff')]
colnames(elnet_all) <- c('from', 'to', 'weight' )
elnet_all <- elnet_all[order(abs(elnet_all$weight), decreasing = TRUE),]
save(elnet_all, file='../Results/Bignet/elnet/elnet_all.obj')
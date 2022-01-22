## Code to generate Gene regulatory networks based on GENIE3, GGM(genet), ARACNE, CLR, Deconvolution and Silencing
source('deps/dircreater.r')
source('deps/wrap_genie3.R')
source('deps/wrap_gennet.r')
source('deps/wrap_minet.r')
source('deps/wrap_decon.r')
source('deps/wrap_silenc.r')
set.seed(2019)
#load replicate wise counts tMM and voom standaradised merch1519arma17trd HTSeq have been filtered for genes with cpm>1 in at least 9 reps
load('../Data/20191112me1915arma17_exp.obj')

#z-scale gene wise
gene_dat <- t(scale(t(me1519arma17_exp)))
#check if transformation has introduced NA values
if(sum(is.na(gene_dat)!=0)) {stop('transformation introduced NA values check expression data')}
#import blast complemented set of TFS from PlnTFDB publication
getmedian <- function(x, gene_dat){ # ADJUSTED FOR REPLICATE INDEX STARTING WITH _X
  #check if data has replicas 
  if (sum(gsub('_[[:digit:]]', '', colnames(gene_dat))==x)>1) {
    return(apply(gene_dat[, gsub('_[[:digit:]]', '', colnames(gene_dat))==x],1, median))
  } else {return(gene_dat[, gsub('_[[:digit:]]', '', colnames(gene_dat))==x])}
} 
gene_dat_median <- sapply(unique(gsub('_[[:digit:]]', '', colnames(gene_dat))), getmedian, gene_dat=gene_dat)


#Import TFs - names should be same as row names of sequencing data
TF <- read.delim('../Data/BLAST_curated_JIN_TF.txt', header =TRUE)
tfs <- TF$final_ID
tfs <- tfs[tfs %in% rownames(gene_dat)]


#Generate GENIE3 network - use median aggregated data
gen3_all <- wrap_genie3(gene_dat_median, tfs=tfs, parallel = 6)
#Generate GGM network - use median aggregated data
genet_all <-wrap_gennet(gene_dat_median, tfs=tfs)
#Generate ARACNE netrowk - use replicate wise data
aracne_all <- wrap_minet(gene_dat, tfs=tfs, algo='ARACNE')
#GENERATE CLR network - use replicate wise data
clr_all <- wrap_minet(gene_dat, tfs=tfs, algo='CLR')
#GEnerate deconvolution network - use median aggregated data
dec_all <- wrap_decon(gene_dat_median,tfs=tfs,silent=TRUE)
#Generate silencing network - use median aggregated data
sil_all <- wrap_silenc(gene_dat_median, tfs=tfs)
#load elastic net results, comment out if elastic net script has not been run 
load('../Results/Bignet/elnet/elnet_all.obj')
#combine results - adapt if elastic net script has not been run
allnet <- list(gen3 = gen3_all[order(abs(gen3_all$weight), decreasing = TRUE),], genet = genet_all[order(abs(genet_all$weight), decreasing = T),],
               elnet = elnet_all[order(abs(elnet_all$weight), decreasing=TRUE),],
               aracne = aracne_all[order(abs(aracne_all$weight), decreasing = TRUE),], 
               clr = clr_all[order(abs(clr_all$weight), decreasing = TRUE),], decon = dec_all[order(abs(dec_all$weight), decreasing=TRUE),],
               silenc = sil_all[order(abs(sil_all$weight), decreasing = TRUE),])#
#final check if there are any G-> XY interactions (non tf as from-node)
print('any non-TF starting nodes')
print(lapply(allnet, function(x) {return(sum(!(x$from %in% tfs)))}))
#Export the unmatched GRNs
dircreater('../Results/Bignet/')
save(allnet, file='../Results/Bignet/allnet.obj')
print('Size of the edge lists returned from the different approaches')
sapply(allnet, dim)
## TF analysis
setwd("~/Chlamy_project_ma/Marius/Program")
source('dircreater.r')
source('Bignet/wrap_elnetv2.r')
source('Bignet/wrap_genie3.R')
source('Bignet/wrap_gennet.r')
source('Bignet/wrap_minet.r')
source('Bignet/wrap_decon.r')
source('Bignet/wrap_silenc.r')
set.seed(2019)
options(stringsAsFactors = FALSE)
#load replicate wise counts tMM and voom standaradised merch1519arma17trd HTSeq have been filtered for genes with cpm>1 in at least 9 reps
load('../Data/20191112me1915arma17_exp.obj')


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

#Import Blast completed TF set of PlnTFDB, JIN(PlantTFDB) and manuall CCM NPQ TFs
TF <- read.delim('../Data/BLAST_curated_JIN_TF.txt', header =TRUE)
tfs <- TF$final_ID
tfs <- tfs[tfs %in% rownames(gene_dat)]
## elastic net run from seperate scritp 'NEtreg_elnet_TF_G_XYZ.R

#use GENIE 3 (decision tree)
gen3_all <- wrap_genie3(gene_dat_median, tfs=tfs, parallel = 3)
genet_all <-wrap_gennet(gene_dat_median, tfs=tfs)
aracne_all <- wrap_minet(gene_dat, tfs=tfs, algo='ARACNE')
clr_all <- wrap_minet(gene_dat, tfs=tfs, algo='CLR')
dec_all <- wrap_decon(gene_dat_median,tfs=tfs,silent=TRUE)
sil_all <- wrap_silenc(gene_dat_median, tfs=tfs)
#load elnetresults
load('../Results/Bignet/me1519am17/TF-G/elnet2/elnet_all.obj')
#combine results
allnet <- list(gen3 = gen3_all[order(abs(gen3_all$weight), decreasing = TRUE),], genet = genet_all[order(abs(genet_all$weight), decreasing = T),],
               elnet = elnet_all[order(abs(elnet_all$weight), decreasing=TRUE),],
               aracne = aracne_all[order(abs(aracne_all$weight), decreasing = TRUE),], 
               clr = clr_all[order(abs(clr_all$weight), decreasing = TRUE),], decon = dec_all[order(abs(dec_all$weight), decreasing=TRUE),],
               silenc = sil_all[order(abs(sil_all$weight), decreasing = TRUE),])#
#final check if there are any G-> XY interactions (non tf as from-node)
print('any non-TF starting nodes')
print(lapply(allnet, function(x) {return(sum(!(x$from %in% tfs)))}))
print(lapply(allnet, head))
dircreater('../Results/Bignet/me1519am17/TF-G')
save(allnet, file='../Results/Bignet/me1519am17/TF-G/allnet.obj')
sapply(allnet, dim)
#use corrmat for silencing and deconvolution, check only edges from TFs
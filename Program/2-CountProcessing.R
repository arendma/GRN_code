# Script combindes read count data from different studies into one gene name based format and 
# then jointly applies TMM and voom normalization of the raw count as implemented in the edgeR package
source('deps/dircreater.r')
source('deps/map55_53.R')
library(edgeR)
library(ggplot2)
library(ggfortify)

###PART 1 Data processing
#load countfiles arma17
load('../Data/arma17_trd.obj')
colnames(arma17) <- paste(rep('am17X', ncol(arma17)), gsub('\\.', '_', colnames(arma17)), sep='')
arma17 <- data.frame(ID5_5=rownames(arma17), arma17)

#load replicate wise counts merch19 HTSeq
load('../Data/merch19.obj')
#'-' is substited for . when converted to data frame therefore substitute '-' for 'neg'
colnames(merch19) <- paste(rep('me19X', ncol(merch19)), gsub('-', 'neg',colnames(merch19)), sep='')
merch19 <- data.frame(ID5_5 =rownames(merch19), merch19)
#load replicate wise counts merch15 HTSeq
merch15 <- read.delim('../Data/GSE71469_ChlamydomonasSynchronousDiurnalExpressionHTSeqCounts.txt', header=TRUE)
#convert 55 to 53 genenames
conv5553 <- map55_53()
merch15_tmp <- merge(conv5553, merch15, by.x='ID5_3', by.y ='Locus.ID..v5.3.1.')
#cleanup and combine in one matrix
merch15_tmp$ID5_3 <- NULL
colnames(merch15_tmp)[2:ncol(merch15_tmp)] <- paste(rep('me15', dim(merch15_tmp)[2]-1), colnames(merch15_tmp)[2:ncol(merch15_tmp)], sep='')
RNA_dat <- merge(arma17, merch15_tmp, by='ID5_5')
RNA_dat <- merge(RNA_dat, merch19, by='ID5_5')
rownames(RNA_dat) <- RNA_dat$ID5_5
RNA_dat <- as.matrix(RNA_dat[2:ncol(RNA_dat)], rownames.force = T)
#select for genes notably expressed in at least 3 conditions e.g. 9 replicates
RNA_fdat <- RNA_dat[rowSums(cpm(RNA_dat)>1) >= 9,]
#clean up cuffdiff condition names
condnames <- factor(gsub('_[[:digit:]]', '', colnames(RNA_fdat)))
#apply standard count data processing
RNalldge <- DGEList(counts=RNA_fdat,group = condnames)
#TMM normalization
RNalldge <- calcNormFactors(RNalldge)
RNalldesign <- model.matrix(~0+condnames)
#dipersion estimate
RNalldge <- estimateDisp(RNalldge, RNalldesign, robust=TRUE)
#voom normalization
RNallv <- voom(RNalldge,RNalldesign,plot=F)
summary(RNallv$E)
## Export the voom normalized (log) read count table so we can work with it later
me1519arma17_exp <- RNallv$E
save(me1519arma17_exp, file='../Data/20191112me1915arma17_exp.obj')

#### PART 2 inital QC analysis 
#generate a light intensity annotation for each sample
samplelight <- rep('LL', 158)
samplelight[unlist(sapply(c('me19Xneg', 'me19X13', 'me15X12', 'me15X13', 'me15X14', 'me15X15', 'me15X16', 'me15X17', 'me15X18', 'me15X19', 'me15X2[[:digit:]]' ), grep, x=colnames(me1519arma17_exp)))] <- 'dark'
samplelight[unlist(sapply(c('am17X5', 'am17X6', 'am17X7', 'am17X8', 'am17X11', 'am17X12', 'am17X13', 'am17X14', 'am17.*HL'), grep, x=colnames(me1519arma17_exp)))] <- 'HL'
samplelight <- factor(samplelight, levels = c('dark', 'LL', 'HL'))
samplelight <- data.frame(samplelight)
rownames(samplelight) <- colnames(me1519arma17_exp)
##Create a PCA of samples
pca1 <- prcomp(t(me1519arma17_exp))
bp <- autoplot(pca1,  label=TRUE, shape=FALSE, label.size=1.5)
dircreater('../Results/Bignet/me1519arma17ana')
ggsave('../Results/Bignet/me1519arma17ana/pcawhole.pdf', bp)

##Create a focussed PCA of our experiments
# design labelst labels for arma17
cultcond <- c(rep(c('HSM', 'HSM', 'HSM', 'Ac', 'Ac', 'Ac'), 7), rep(c('HSM_phot', 'HSM', 'HSM_phot', 'HSM'), each=3))
lightcond <- c(rep('LL', 6), rep(c(rep('LL', 6),rep('HL',12)), 2), rep('HL',6), rep('LL',6))
timep <- c(rep(c(15, 240, 255,300,1440, 1455, 1500), each=6), rep(-1, 12))
timelight <- factor(c(rep(c('LL', 'LL', 'HL_255', 'HL_300', 'LL', 'HL_1455', 'HL_1500'), each=6), rep(c('HL_phot', 'HL_WT', 'LL_phot', 'LL_WT'), each=3)), levels=c('LL', 'HL_255', 'HL_300', 'HL_1455', 'HL_1500', 'HL_phot', 'HL_WT', 'LL_phot', 'LL_WT'))
lab <- gsub('_0', '',paste(cultcond, timep, sep='_'))
arma17cond <- data.frame(cultcond, lightcond, timep, timelight, lab)
rownames(arma17cond) <- grep('am17', colnames(me1519arma17_exp), value=TRUE)
arma17_exp <- me1519arma17_exp[,grep('am17', colnames(me1519arma17_exp))]
pca2 <- prcomp(t(arma17_exp))
bp2 <- autoplot(pca2,  data=arma17cond, label.label='lab' , shape='cultcond', colour='lightcond' ,size=7) + theme_bw()+theme(text=element_text(size=25)) + labs(colour='Light treatment', shape='Algae culture')
ggsave('../Results/Bignet/me1519arma17ana/pcaarma.pdf', bp2, height = 8, width=12, useDingbats=FALSE)

#recuce pca coordinates to only show timeline
pca3 <- pca2
pca3$x <- pca3$x[1:42,]
ma17cond<- arma17cond[1:42,]
bp3 <- autoplot(pca3, data= ma17cond, shape='cultcond', colour='timelight', size=7) + theme_bw() +theme(text=element_text(size=25))+ labs(colour='Light treatment + time', shape='Algae culture')
ggsave('../Results/Bignet/me1519arma17ana/pcama17.pdf', bp3, height=8, width=12 ,useDingbats=FALSE)

##Create a heat map of sample correlations
library(amap)
library(pheatmap)
corm=cor(me1519arma17_exp)
d = Dist(x=t(me1519arma17_exp), method = 'correlation')
custome_color <- list(samplelight = c(dark = 'black', LL= 'grey50', HL='grey90'))
pdf('../Results/Bignet/me1519arma17ana/realcorrsamples2.pdf', width=30, height=30)
pheatmap(as.matrix(corm), clustering_distance_rows = d, clustering_distance_cols = d, annotation_row = samplelight, annotation_col = samplelight, annotation_colors=custome_color)
dev.off()
pdf('../Results/Bignet/me1519arma17ana/realcorrsamples_small.pdf', width=8, height=8)
pheatmap(as.matrix(corm), clustering_distance_rows = d, clustering_distance_cols = d, annotation_row = samplelight, annotation_col = samplelight, annotation_colors=custome_color, show_rownames=FALSE, show_colnames=FALSE)
dev.off()

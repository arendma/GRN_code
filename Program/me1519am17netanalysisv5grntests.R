#Script to build consenesus network with all networks present and consensus networ w/o ARACNE/silencing
setwd("C:/Users/Suiram/Nextcloud/PhD/PHOT_Project/Chlamy_project_ma/Marius/Program/")
source("comp_edgev3.r")
source('dircreater.r')
source('map55_4.R')
source('Bignet/netwk_compare.R')
source('Bignet/networkmerge.R')
library(igraph)
conversion <- map55_4()
#looad consensus network 
load('../Results/Bignet/me1519am17/TF-G/allnet.obj')
me1519am17allnet <- allnet
rm(allnet)
#set limit of used genes to 10% of all possible edges
topx <- round(0.1*nrow(me1519am17allnet[['decon']]))

#adapt pulicable names
names(me1519am17allnet) = c('GENIE3', 'GGM', 'elastic_net', 'ARACNE', 'CLR', 'Deconvolution', 'Silencing')
#keep ARACNE and Silencing for comparison
#do a consensus network containing all methods
me1519am17allconsensus <- networkmerge(me1519am17allnet, topx=topx, resdir='../Results/Bignet/me1519am17/TF-G/allnet/')



me1519am17final=me1519am17allnet
me1519am17final[c('ARACNE', 'Silencing')]=NULL
### Build the Consensus networks
me1519am17consensus <-networkmerge(me1519am17final, topx=topx, resdir='../Results/Bignet/me1519am17/TF-G/10pcfinal/')

#old code
if (FALSE) {
truinteract <- read.delim('../Data/regulation_from_motif_Cre.txt', header = FALSE)
truinteract <- data.frame(from = truinteract[,1], to = truinteract[,3])
#save in as consensusfile without weight to be able to run netwk_comp
write.table(truinteract, file='../Data/regulation_from_motif_Cre_netwk.tab', row.names=FALSE, col.names=TRUE, sep='\t')
print(paste('Of the ', length(unique(truinteract[,'from'])), 'TFs present in the JIN network', length(intersect(me1519am17consensus[,'from'], truinteract[,'from'])), 
          'TFs are present in the network', collapse = ''))
#calculate number of found interactions
found_interacts <- sum(comp_edge(me1519am17consensus[,c('from', 'to')], truinteract))
print(paste('Of the ', nrow(truinteract), 'interactions in the JIN network' , found_interacts,  '(', 
            (found_interacts/nrow(truinteract))*100, '%) are found in the consensus network', collapse=''))
#reduce reported edges to be able to compare methods
for(approach in names(me1519am17allnet)[names(me1519am17allnet)!='ARACNE']) {
  me1519am17allnet[[approach]] <- me1519am17allnet[[approach]][1:topx,]
}
netwk_compare('../Results/Bignet/me1519am17/TF-G/allnet/consensus.tab', netwks = me1519am17allnet,topx=100000, resdir = '../Results/Bignet/me1519am17/TF-G/allnet/allconsensus/')
netwk_compare('../Results/Bignet/me1519am17/TF-G/allnet/consensus.tab', netwks = me1519am17allnet,topx=500000, resdir = '../Results/Bignet/me1519am17/TF-G/allnet/allconsensus5e7/')
netwk_compare('../Data/regulation_from_motif_Cre_netwk.tab', netwks = me1519am17allnet,topx=nrow(truinteract), resdir = '../Results/Bignet/me1519am17/TF-G/allnet/motifgrn/')

#import CCM regulatory network adjacency matrix from winck 2013 publication
CCM_adjmat <- read.delim('../Data/CCM_GRN_adjmat_winck2013.txt', row.names=1)
CCM_adjmat <- as.matrix(CCM_adjmat)
colnames(CCM_adjmat) <- gsub('X', '', colnames(CCM_adjmat))
#try to rescue JGI v4.0 annotations
newcoln <- unlist(sapply(colnames(CCM_adjmat), match_4_55, convetab=conversion))
newrown <- unlist(sapply(rownames(CCM_adjmat), match_4_55, convetab=conversion))

nCCM_adjmat<- CCM_adjmat[which(!(is.na(newcoln))),which(!(is.na(newrown)))]
colnames(nCCM_adjmat) <- newcoln[!(is.na(newcoln))]
rownames(nCCM_adjmat) <- newrown[!(is.na(newrown))]

CCM_g <- graph_from_adjacency_matrix(nCCM_adjmat, weighted=TRUE)

#get edgelist from graph
CCM_g_elist <- as_edgelist(CCM_g)
colnames(CCM_g_elist) <- c('from', 'to')

print(paste('Of the ', ncol(nCCM_adjmat), 'TFs in the CCM TF-TF network' ,length(intersect(CCM_g_elist[,'from'], me1519am17consensus[,'from'])),
            'are present in the consensus network.', collapse=''))
found_interacts2 <- sum(comp_edge(CCM_g_elist, (me1519am17consensus[, c('from', 'to')])))

print(paste('Of the ', nrow(CCM_g_elist), 'interactions in the CCM TF-TF network ' , found_interacts2, '(', 
            (found_interacts2/nrow(CCM_g_elist))*100, '%) are present in the consensus network', collapse=''))
consensusnet <- graph_from_data_frame(me1519am17consensus[, c('from', 'to')], directed=T)

compnet <- difference(CCM_g, consensusnet)

E###################################################################################################
#Import network from elastic regression
load('../Results/Save20200108/Transcriptome/TF_network/TF-G/elnet_all.obj')
ph_elnet <- elnet_all
rm(elnet_all)
#Import network from GENIE3
load('../Results/Save20200108/Transcriptome/TF_network/TF-G/gen3_all.obj')
ph_gen3 <- gen3_all
#Since gen3 returns all possible edges only take 1top 10% for photmot network and for bignet the same number of edges as contained in the consensus network
ph_gen3top <- ph_gen3[1:round(nrow(ph_gen3)*0.1),]
rm(gen3_all)

me1519am17allconsensus <- networkmerge(me1519am17allnet, topx=topx, resdir='../Results/Bignet/me1519am17/TF-G/10pcfinal/')



allnet <- c(me1519am17allnet, consensus=list(me1519am17consensus[,c('from', 'to', 'meanrank')]),ph_elnet=list(ph_elnet), ph_gen3=list(ph_gen3))

litedges <- read.delim('../Data/TF-Gnpqccm.txt')
helpf1 <- function(x) {
  #takes an egdge as a character vectro of length 2 and looks for it in all rankings in me1519am17allnet. If it is not present it returns NA
  res <- numeric(0)
  for (method in names(allnet)){
    rank <- which(comp_edge(allnet[[method]][,1:2],data.frame(x[1], x[2])))
    if (length(rank)==0) {res <- c(res, NA)}
    else if(length(rank)==1) {res <- c(res, rank)}
    else if(length(rank)>1) {stop('duplicated edges detected :(')}
  }
  names(res) <- names(allnet)
  return(res)
}
method_ranks <- t(apply(litedges[,c('TFID5_5', "targetID5_5")], 1, helpf1))
#gen3_ranks <- apply(npqedges, 1, function(x) {return(which(comp_edge(me1519am17allnet[['gen3']][,1:2],data.frame(x[1], x[2]))))})
#decon_ranks <-apply(npqedges, 1, function(x) {return(which(comp_edge(me1519am17allnet[['decon']][,1:2],data.frame(x[1], x[2]))))})
npqedgesres <- data.frame(litedges[,c('TFname', "targetname")], method_ranks)
#normalize to edgenumber reported by method
if(identical(names(allnet),names(sapply(allnet, dim)[1,]))) {
  npqedgesresrel <- npqedgesres
  npqedgesresrel[,3:ncol(npqedgesresrel)]<- t(t(npqedgesres[,3:ncol(npqedgesresrel)])/sapply(allnet, dim)[1,])
} else {stop('ERROR!order of methods in literature interaction dataframe and network object are not identical.')}
library(pheatmap)
pdf('../Results/Bignet/me1519am17/TF-G/10pcfinal/litinteraction_aracnesil.pdf')
pheatmap(as.matrix(npqedgesresrel[,3:ncol(npqedgesresrel)]), cluster_rows=FALSE, cluster_cols=FALSE, fontsize=14, labels_row = paste(npqedgesres[,1],npqedgesres[,2], sep='->'), display_numbers=npqedgesres[,3:ncol(npqedgesresrel)],fontsize_number=7)
dev.off()

}

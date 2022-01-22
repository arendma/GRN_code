#This code loads the previously created GRNs and matches the edges to obtain a consensus network
source("deps/comp_edgev3.r")
source('deps/dircreater.r')
source('deps/networkmerge.R')
library(igraph)
#conversion <- map55_4()
#looad consensus network 
load('../Results/Bignet/allnet.obj')
me1519am17allnet <- allnet
rm(allnet)
#set limit of used genes to 10% of all possible edges
topx <- round(0.1*nrow(me1519am17allnet[['decon']]))
#adapt pulicable names
names(me1519am17allnet) = c('GENIE3', 'GGM', 'elastic_net', 'ARACNE', 'CLR', 'Deconvolution', 'Silencing')
#keep ARACNE and Silencing for comparison
#do a consensus network containing all methods
me1519am17allconsensus <- networkmerge(me1519am17allnet, topx=topx, resdir='../Results/Bignet/allnet/')


#Create the final consensus network omitting ARACNE and Silencing results
me1519am17final=me1519am17allnet
me1519am17final[c('ARACNE', 'Silencing')]=NULL
### Build the Consensus networks
me1519am17consensus <-networkmerge(me1519am17final, topx=topx, resdir='../Results/Bignet/10pcfinal/')
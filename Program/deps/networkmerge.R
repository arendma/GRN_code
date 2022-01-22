networkmerge <- function(netwks, topx, resdir) {
  source('dircreater.r')
  source('comp_edge.r')
  #given a list of adjacency datframes(rows=edges, columns: ['from']=from node, ['to']=to node,  ['weight']=edge weight)
  #and merges the topx notes
  dircreater(resdir)
  checkdim <- sapply(netwks, dim)
  if (max(checkdim) < topx) {stop('ERROR: topx > max(lapply(netwks, dim)[1,] ')}
  helpf1 <- function(df){
    if (any(df$weight!=df$weight[order(abs(df$weight), decreasing = TRUE)])){
      print(paste('Warning unsorted edgelists detected... sorting edges according to weight'))
      return(data.frame(df[order(abs(df$weight),decreasing = TRUE)], ranks=1:dim(df1)[1]))
    } else {
      return(data.frame(df, ranks=1:dim(df)[1]))
    }
  }
  netwks <- lapply(netwks, helpf1)
  #to reduce comparisons nescessary reorder model list, so that model with maximum amount of edges ist first model
  
  netwks <- netwks[order(checkdim[1,], decreasing=TRUE)]
  
  #aggregate list of all edges
  consensus <- netwks[[1]][1:topx,c('from', 'to', 'ranks')]
  colnames(consensus) <- c('from', 'to', names(netwks)[1])
  for (i in 2:length(netwks)) {
    mdl <- netwks[[i]]
    if (dim(mdl)[1]<topx) {
      limit <- dim(mdl)[1]
    } else {
      limit <- topx
    }
    consensus <- merge(consensus, mdl[1:limit, c('from', 'to', 'ranks')], by=c('from', 'to'), all=TRUE)
    colnames(consensus)[dim(consensus)[2]] <- names(netwks)[i]
  }
  
  #If edge is not found by model give rank topx+1
  consensus[,3:dim(consensus)[2]] <- apply(consensus[,3:dim(consensus)[2]], 2, function(x){
    x[is.na(x)]<- topx+1
    return(x)
    })
  #generate mean ranks
  consensus <-data.frame(consensus, meanrank = apply(consensus[,3:dim(consensus)[2]],1, mean))
  consensus <- consensus[order(consensus$meanrank),]
  rownames(consensus) <- apply(consensus[,1:2], 1, paste, collapse=' ')
  write.table(consensus, file=file.path(resdir,'consensus.tab'), sep='\t',)
  return(consensus)
  if (FALSE) {#old code
    ##ATTENTION THIS CODE PRODUCED DUPLICATES....
    if(FALSE){
      for(l in 1:length(netwks)) {
        if (l==1) {
          #intialis consensus dataframe
          consensus <- cbind(netwks[[l]][1:topx,1:2],1:topx, matrix(ncol=length(netwks)-1, nrow=topx))
          colnames(consensus)<- c('from', 'to', names(netwks))
        } else {
          mdl <- netwks[[l]][1:topx,]
          for (i in 1:dim(mdl)[1]) {
            test_e <- mdl[i,]
            if (test_e[1] %in% consensus[,'from']) {
              #report save rank if edge is contained
              if (test_e[2] %in% consensus[consensus$from %in% test_e[1],'to']) {
                consensus[consensus[,'from'] %in% test_e[1] & consensus[,'to'] %in% test_e[2], l+2] <- i 
              }
              else {
                #build new entry is not present in consensus data frame
                new_e <- c(test_e[1:2], rep(NA, length(netwks)))
                names(new_e) <- colnames(consensus)
                new_e[l+2] <- i
                #append
                consensus <- rbind(consensus, new_e, make.row.names = FALSE)
              }
              
            } else{
              #same as upper else loop
              ew_e <- c(test_e[1:2], rep(NA, length(netwks)))
              names(new_e) <- colnames(consensus)
              new_e[l+2] <- i
              consensus <- rbind(consensus, new_e, make.row.names=FALSE)
            }
          }
        }
      }
    } 
  #parse data from http://plantregmap.cbi.pku.edu.cn/download.php
  truinteract <- read.delim('../Data/regulation_from_motif_Cre.txt', header = FALSE)
  truinteract <- data.frame(from = truinteract[,1], to = truinteract[,3])
  # select only tfs found in our data
  print(paste('Number of motif regulation TFS from jin PlantTF found in our data:', sum(unique(truinteract$from) %in% consensus$from)))
  #only keep TFis included in our set
  truinteract <- truinteract[truinteract$from %in% consensus$from,]

  print(paste('rel. amount of motif predicted interactions for these found:', sum(comp_edge(truinteract, consensus[,1:2]))))
  sum(comp_edge(truinteract, consensus[,1:2]))
  #test function
  test<- truinteract[1:100,]
  #check for found interactions
  print(paste('interactions found in motif scanning based data', sum(comp_edge(rbind(consensus_e[,1:2]),truinteract))))
  #compare consensus with single methods
  metcomp <- matrix(nrow=10*length(netwks), ncol=3)
  for (i in 1:length(netwks)) {
    for (j in 1:10) {
      subs <- netwks[[i]][1:round(dim(netwks[[i]])[1]*(j/10)),1:2]
      metcomp[(i*10-10+j),1] <- sum(comp_edge(subs,consensus))
      metcomp[(i*10-10+j),2] <- j/10
      metcomp[(i*10-10+j),3] <- names(allnet)[i]
    }
  }
  colnames(metcomp) <- c('found', 'fraction', 'method')
  metcomp <- data.frame(found = as.numeric(metcomp[,1]), fraction = as.numeric(metcomp[,2]), method =metcomp[,3])
  
  ggplot(data=metcomp, aes(x=fraction, y=found, group=method, colour=method)) + 
    geom_line(size=1.5) + xlab('% of found interactions') + ylab('number of consensus interactions')
  geom_point(size=3, fill="white") 
  ggsave(file.path(resdir, 'internal_consensus_comp.pdf'))
  }
}
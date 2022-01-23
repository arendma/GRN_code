networkmerge <- function(netwks, topx, resdir) {
  source('deps/dircreater.r')
  source('deps/comp_edgev3.r')
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
}
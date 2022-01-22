map55_53 <- function() {
  #function returns a unambigous dataframe with to columns to convert from 5.3 to 5.5 annotation
  prtnclean <- function(prtn) {
    return(gsub("\\.t.*", "", prtn))
  }
  con <- file('../Data/ChlamydomonasTranscriptNameConversionBetweenReleases.Mch12b.txt', open = 'r')
  gn <- list()
  while (TRUE) {
    #parser for conversion  file with 25 char spacing for each colum
    line = readLines(con, n=1)
    if (length(line)==0) {break}
    nchars <- nchar(line)
    if (nchars <50) {next}
    start <- seq(1,nchars-24, 25 )
    end <- seq(24, nchars, 25)
    ldat <- substring(line,start, end)
    ldat <- gsub(' ', '' , ldat)
    gn <- rbind(gn, gsub('--', NA, ldat))
  }
  close(con)  
  #assign headers
  colnames(gn) <- gn[1,]
  gn <- gn[-1,]
  conversion <- data.frame(ID5_5=prtnclean(gn[,'#5.5']), ID5_3 = prtnclean(gn[,'5.3.1']))
  
  collapsegene <- function(gndf) {
    resgndf <- gndf
    #take th data frame 
    #for each 5.5 gene id that appears more than once
    for (dupgn in unique(gndf[duplicated(gndf[,1]),1])) {
      idx <- which(resgndf[,1]==dupgn)
      #print(resgndf[idx,])
      #check if all locusses are linked to the same gene name
      if(length(unique(resgndf[idx,2]))>1) {
        #if not concatenate the v.4.0 annotions to one string and link it to the first occurence of the locus
        print('yes')
        resgndf[idx[1],2] <- paste(as.character(unique(resgndf[idx,2])), collapse=':')
      }
      #remove all other locus occurences
      resgndf <- resgndf[-(idx[2:length(idx)]),]
    }
    return(resgndf)
  }
  conversion <- collapsegene(conversion)
  n_duplicated <- dim(conversion)[1]
  #remove entries of duplicated 5.3 annotations that map to different 5.5 entries
  for (dupgene in unique(conversion$ID5_3[duplicated(conversion$ID5_3)])) {
    conversion <- conversion[-(which(conversion$ID5_3==dupgene)),]
  }
  print(paste(as.character(n_duplicated-dim(conversion)[1]), ' genes removed due to ambigous asignment', sep=''))
  return(conversion)
}
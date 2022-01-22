wrap_elnet <- function (dat, resdir, thrsh, tfs=FALSE, parallel=FALSE, ...) {
  #takes a matrix of gene expression values and calculates elastic net regression coefficients for all gene combinations
  #gene names must be specified in rownames(data)
  source('deps/TFelnetv6.r')
  source('deps/dircreater.r')
  library(parallel)
  dircreater(resdir)
  reg_res <- data.frame()
  #predict all genes
  if(tfs==FALSE) {
  #serial
  if (parallel==1) {
    for (i in 1:length(rownames(dat))) {
      reg_res <- rbind(reg_res, TFelnet(y=as.numeric(dat[i,]), x=t(dat[-(i),]), goi_n=rownames(dat)[i], resdir=resdir, ...))
    }
    new_res <- reg_res
  }
  #parallel
  if (parallel>1) {
    ncl <- detectCores()
    if (parallel>ncl) {stop('more threads than cores available!!')}
    cl <- makeCluster(parallel)
    clusterEvalQ(cl, {
      source('deps/TFelnetv6.r')
    })
    f <- function(i) {
      return(TFelnet(y=as.numeric(dat[i,]), x=t(dat[-(i),]), goi_n=rownames(dat)[i], resdir=resdir, ...))
    }
    #normal lapply -- works
    #reg_res <- lapply(X=1:dim(dat)[1],FUN=f)
    #parallel
    reg_res <- parLapply(cl,X=1:dim(dat)[1],fun=f)
    #close clusters
    stopCluster(cl)
    for (i in 2:length(reg_res)) {
      if (i==2) {
        new_res <- rbind(reg_res[[1]], reg_res[[i]], make.row.names=FALSE)
      }
      else {
        new_res <- rbind(new_res, reg_res[[i]], make.row.names=FALSE)
      }
    }
  }
  }
  #With tf list as predictors
  else {
    if (parallel==1) {
      predictor<- dat[rownames(dat) %in% tfs,]
    for (i in 1:length(rownames(dat))) {
      if (rownames(dat)[i] %in% rownames(predictor)) {
      xs <- predictor[-which(rownames(predictor)==rownames(dat)[i]),]
      }
      else {
        xs <- predictor
      }
      reg_res <- rbind(reg_res, TFelnet(y=as.numeric(dat[i,]), x=t(xs), goi_n=rownames(dat)[i], resdir=resdir, ...))
    }
    new_res <- reg_res
  }
    #parallel
    if (parallel>1) {
      ncl <- detectCores()
      if (parallel>ncl) {stop('more threads than cores available!!')}
      cl <- makeCluster(parallel)
      clusterEvalQ(cl, {
        source('deps/TFelnetv6.r')
      })
      f <- function(i,tfs) {
        predictor<- dat[rownames(dat) %in% tfs,]
        if (rownames(dat)[i] %in% rownames(predictor)) {
          xs <- predictor[-which(rownames(predictor)==rownames(dat)[i]),]
        }
        else{
          xs <- predictor
        }
        return(TFelnet(y=as.numeric(dat[i,]), x=t(xs), goi_n=rownames(dat)[i], resdir=resdir, ...))
      }
      #normal lapply -- works
      #reg_res <- lapply(X=1:dim(dat)[1],FUN=f)
      #parallel
      reg_res <- parLapply(cl,X=1:dim(dat)[1],fun=f, tfs)
      #close clusters
      stopCluster(cl)
      for (i in 2:length(reg_res)) {
        if (i==2) {
          new_res <- rbind(reg_res[[1]], reg_res[[i]], make.row.names=FALSE)
        }
        else {
          new_res <- rbind(new_res, reg_res[[i]], make.row.names=FALSE)
        }
      }
    }
    }
  #Filter bad models
  relvarsort <- sort(unique(new_res$relvar))
  lost <- sum(relvarsort<thrsh)/length(relvarsort)
  print(paste('Mit Threshold R2>', as.character(thrsh),  'wurden', as.character(lost), '% der modelle verworfen'))
  new_res <- new_res[new_res$relvar>thrsh,]
  #Filter non interacting genes
  new_res <- new_res[abs(new_res$coefficients)>0,]
  #globally weight
  new_res <- data.frame(new_res, rel.coeff = new_res$coefficients/max(abs(new_res$coefficients)))
  reg_res <- new_res
  save(reg_res, file=file.path(resdir, 'elnet.obj'))
  return(reg_res)
}
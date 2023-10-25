TFelnet <- function(y, x, goi_n, resdir, intercept=FALSE, K=6, report=FALSE) {
  set.seed(2019)
  #do elastic net regression for gene of interest profile
  #find optimal lambda and s by crossvalidation and return coefficients for that conditions
  #if data is not centralized intercept should be set two TRUE
  source('deps/genenamecleanupv2.r')
  source('deps/dircreater.r')
  library(elasticnet)
  library(ggplot2)
  #import gene names
  phyto_gn <- genenamecleanup()
  #look for gene of interest in gene names - removed because genen names are not unique
  #goi_n <- gn_match(goi_n, phyto_gn)
  if(report){
    if (!(substr(goi_n, 1,3)=='Cre') & dir.exists(file.path(resdir, goi_n))) {
      print('Warning! filename does already exist, maybe check for doubled genen names')
    }
    dircreater(file.path(resdir, goi_n))
  }
  #center x values rowwise (per condition)
  #x <- t(scale(t(x), center=TRUE, scale=FALSE))
  # set lambda [lambda2 the quadratic punishment]
  lambdas <- c(0,0.001,.01,.05,.1,.5,1,1.5,2,10,100)
  #innitalise matrix to store values for fit with lowest prediction error for each lambda
  cv_res <- matrix(nrow=4, ncol=length(lambdas))
  rownames(cv_res) <- c('lambda', 's', 'cv', 'cv_er')
  cv_res[1,] <- lambdas
  #do the crossvalidation
  for (i in seq(1, length(lambdas))) {
    cv_lm <- cv.enet(x=x, K=K, y=y, s=seq(0.1,1, 0.1), mode='fraction', lambda=lambdas[i], intercept=intercept)
    idx_min <- which.min(cv_lm$cv)
    cv_res[2:4,i] <- c(cv_lm$s[idx_min], cv_lm$cv[idx_min], cv_lm$cv.error[idx_min])
  }
  #only keep res with <0.91
  cv_res <- cv_res[, cv_res['s',]<1 & cv_res['s',] >0.09]
  #plot cv vs lambda
  if(report) {
    p<- ggplot(data.frame(t(cv_res)), aes(x=lambda, y=cv)) + 
      geom_point() +  scale_x_log10() + scale_y_log10() +
      geom_errorbar(aes(ymin=cv-cv_er, ymax=cv+cv_er), width=.2)
    ggsave(file.path(resdir, goi_n, 'cvmin_lamb_s.pdf'))
  }
  #select lambda with lowest cv and build model
  lm <- enet(x=x, y=y, lambda=cv_res[1, which.min(cv_res['cv',])], intercept=intercept)
  if(report) {
  save(lm, file=file.path(resdir, goi_n, 'lm.obj'))
  }
  if(report) {
    pdf(file.path(resdir, goi_n, 'betavs_s.pdf'))
    plot(lm)
    abline(v=cv_res['s', which.min(cv_res['cv',])], col='red', lty='dashed')
    dev.off()
  }
  #select s showing lowest cv and generate coefficients
  beta_lm <- predict(lm, s=cv_res['s', which.min(cv_res['cv',])], type='coefficients', mode ='fraction')
  x_lm <- predict(lm, newx= x, s=cv_res['s', which.min(cv_res['cv',])], type='fit', mode ='fraction')
  regulr2 <- function(obs, fit) {
    return(1-(sum((obs-fit)^2)/sum((obs-mean(obs))^2)))
  }
  relvar <- regulr2(y,x_lm$fit)
  beta_lm_df <- data.frame(predicted = rep(goi_n, length(beta_lm$coefficients)), Gene.ID = names(beta_lm$coefficients), coefficients=beta_lm$coefficients, 
                           relvar=rep(relvar, length(beta_lm$coefficients)), 
                           s=rep(cv_res['s', which.min(cv_res['cv',])], length(beta_lm$coefficients)), lambda=rep(cv_res['lambda', which.min(cv_res['cv',])], length(beta_lm$coefficients)),
                           cv=rep(cv_res['cv', which.min(cv_res['cv',])], length(beta_lm$coefficients)))
  beta_lm_df <- beta_lm_df[order(abs(beta_lm_df$coefficients), decreasing=TRUE),]
  #write.table(beta_lm_df, file = file.path(resdir, goi_n, 'coeff.txt'))
  #p2 <- ggplot(beta_lm_df[abs(beta_lm_df$coefficients) >0.05,], aes(x=Gene.Name, y=coefficients)) + geom_bar(stat='identity')+
  #  theme(axis.text.x = element_text(angle = 45))
  #ggsave(file.path(resdir, goi_n, 'coeff.pdf'))
  return(beta_lm_df)
} 
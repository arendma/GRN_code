Deconvolution <- function(G, silent)
{
  print(isSymmetric(G))
  
  print(all(diag(G)==0))
  print(all(diag(G)==1))
  
  us_G_obs <- G
  
  Beta=0.9
  
  repeat{
    # 1- linear Scaling Step
    G_obs <- linear_scaling(us_G_obs, Beta=Beta) 
    
    # 2- Decomposition Step
    
    r <- try(eigen(G_obs) ,silent = F)
    U <- r$vectors
    lam_obs <- r$values
    
    # 3- Deconvolution Step
    
    lam_dir <- lam_obs/(lam_obs+1)
    
    Lmbd_dir <- diag(lam_dir)
    print(max(abs(Lmbd_dir)))
    if(max(abs(Lmbd_dir))>1)
    {
      if (silent){
        Beta <- Beta-0.05
        print(paste('Beta =', Beta))
      }
      else {
      print(paste("please give a beta less than",Beta,sep=" "))
      Beta <- as.numeric(readline())
      print(Beta)
      }
    }
    else break
  }
  
  sU <- solve(U) 
  
  G_dir <- U %*% Lmbd_dir %*% sU
  colnames(G_dir) <- colnames(us_G_obs)
  rownames(G_dir) <- rownames(us_G_obs)
  
  diag(G_dir) <- 0
  
  # reRANK
  #G_dir <- abs(G_dir)/max(abs(G_dir))
  
  return(G_dir)
}

linear_scaling <- function(us_G_obs, Beta )
{
  G1 <- us_G_obs
  
  r_G1 <- try(eigen(G1) ,silent = F)
  
  U_us <- r_G1$vectors; 
  lam_obs_us <- r_G1$values
  
  lam_obs_us_pos <- max(lam_obs_us)
  lam_obs_us_neg <- min(lam_obs_us)
  
  alpha = max(Beta/((1-Beta)*lam_obs_us_pos),-1*Beta/((1+Beta)*lam_obs_us_neg))
  
  return(alpha * us_G_obs)
  
}
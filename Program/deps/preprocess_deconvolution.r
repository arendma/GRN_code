Dpreprocess <- function(G, tfs, nontfs, silent)
{
  
  N <- dim(G)[2]
  
  # correlation matrix
  cor_matd <- G
#   cor_matd[nontfs,] = 0
  
  G_obs <- cor_matd 
  diag(G_obs) <- 1
  
  G_dir <- Deconvolution(G_obs, silent)
#   G_dir[nontfs,] = 0
  diag(G_dir) <- 0
  
  return(G_dir)
  
}


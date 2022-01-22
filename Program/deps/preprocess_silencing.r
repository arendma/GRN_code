# Besme Allahe Rahmane Rahim
# Khodaya be omide to na be omide khalghe roozegar      Ya HAGH


Spreprocess <- function(G, tfs, nontfs)
{
  
  N <- dim(G)[2]
  
  # correlation matrix
  cor_mat <- G
  
#   cor_mat <- rbind(cor_mat[tfs,],cor_mat[nontfs,])
  cor_mat[nontfs,] = 0
  
#   cor_mat <- cbind(cor_mat[,tfs],cor_mat[,nontfs])
  
  print(all(cor_mat[nontfs,]==0))
  
  diag(cor_mat) <- 1
  
  S <- Silencing(cor_mat)

  return(S)
  
}

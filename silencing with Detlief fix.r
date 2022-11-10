Silencing <- function(G)
{
  diag(G)=1
  
  I <- diag(ncol(G))
  
  D <- diag((G - I) %*% G)
  
  c_inv <- try(solve(G),silent = F)
  print(class(c_inv))
  if (class(c_inv)[1] == "try-error")
    c_inv <- mpinv(G)
  
  S <- (G - I + diag(D)) %*% c_inv#solve(G)
  
  diag(S) <- 0
  
  #S = abs(S)/max(abs(S))
  
  return(S)
}
mpinv <- function(A, eps = 1e-13) {
  s <- svd(A)
  e <- s$d
  e[e > eps] <- 1/e[e > eps]
  return(s$v %*% diag(e) %*% t(s$u))
}
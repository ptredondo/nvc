kencoh.eigen <- function(Ryx, Ryy.inv, Rxx.inv, vectors = TRUE) {
  
  Rxy <- t(Ryx)
  
  ## Matrix determining canonical coherence
  A <- Ryy.inv %*% Ryx %*% Rxx.inv %*% Rxy
  
  ## Only eigenvalues during lag search
  if (!vectors) {
    eg <- eigen(A, only.values = TRUE)$values
    return(eg[1])
  }
  
  ## Full decomposition at the selected lag
  eg.y <- eigen(A)
  A2 <- Rxx.inv %*% Rxy %*% Ryy.inv %*% Ryx
  eg.x <- eigen(A2)
  
  return(list(
    eig.y = eg.y$values, eigvec.y = eg.y$vectors,
    eig.x = eg.x$values, eigvec.x = eg.x$vectors,
    A = A
  ))
}


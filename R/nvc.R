nvc <- function(x, y, perm.inv = TRUE, nperm = NULL,
                blength = 100, freq = NULL, b.scale = FALSE,
                nboot = 0, ncores = 1) {
  
  Tlength <- nrow(y)
  bnum <- Tlength / blength
  
  ## Block periodograms for X
  x.f <- array(NA, dim = c(bnum, ncol(x), blength / 2))
  
  for (p in 1:ncol(x)) {
    for (b in 1:bnum) {
      if (b.scale) {
        x.f[b, p, ] <- LSTS::periodogram(
          scale(x[blength * (b - 1) + 1:blength, p]),
          plot = FALSE
        )[[1]]
      } else {
        x.f[b, p, ] <- LSTS::periodogram(
          x[blength * (b - 1) + 1:blength, p],
          plot = FALSE
        )[[1]]
      }
    }
  }
  
  ## Block periodograms for Y
  y.f <- array(NA, dim = c(bnum, ncol(y), blength / 2))
  
  for (q in 1:ncol(y)) {
    for (b in 1:bnum) {
      if (b.scale) {
        y.f[b, q, ] <- LSTS::periodogram(
          scale(y[blength * (b - 1) + 1:blength, q]),
          plot = FALSE
        )[[1]]
      } else {
        y.f[b, q, ] <- LSTS::periodogram(
          y[blength * (b - 1) + 1:blength, q],
          plot = FALSE
        )[[1]]
      }
    }
  }
  
  ## Fourier frequencies
  if (is.null(freq)) {
    freq <- seq(0, 0.5, length = blength / 2 + 1)[-1]
  }
  
  ## Initialize NVC estimates
  nvc.est <- rep(NA, length(freq))
  
  ## NVC estimation with permutation testing
  if (nboot > 0) {
    
    pval <- rep(NA, length(freq))
    
    for (f in seq_along(freq)) {
      res <- xi.vec(
        x = x.f[, , f],
        y = y.f[, , f],
        perm.inv = perm.inv,
        nperm = nperm,
        nboot = nboot,
        ncores = ncores
      )
      
      nvc.est[f] <- res$xi.vec
      pval[f] <- res$pval
    }
    
    return(list(
      freq = freq,
      nvc = nvc.est,
      pval = pval,
      pval.adj = stats::p.adjust(pval, method = "BH")
    ))
    
  } else {
    
    ## NVC estimation without permutation testing
    for (f in seq_along(freq)) {
      nvc.est[f] <- xi.vec(
        x = x.f[, , f],
        y = y.f[, , f],
        perm.inv = perm.inv,
        nperm = nperm,
        nboot = nboot,
        ncores = ncores
      )
    }
    
    return(list(
      freq = freq,
      nvc = nvc.est
    ))
  }
}




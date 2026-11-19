cancoh <- function(x, y, samprate = 1, blength = 100,
                   demean = TRUE, ridge = 0, tol = 1e-10) {
  
  ## Input checks
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  if (nrow(x) != nrow(y)) {
    stop("x and y must have the same number of observations.")
  }
  
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(y)
  
  if (blength > n) {
    stop("blength cannot exceed the number of observations.")
  }
  
  B <- blength
  starts <- seq(from = 1, to = n - B + 1, by = B)
  bnum <- length(starts)
  
  if (bnum < 1) {
    stop("No complete blocks are available.")
  }
  
  ## Combine X and Y
  Z <- cbind(x, y)
  d <- p + q
  
  ## Hann window and normalization
  j <- 0:(B - 1)
  w <- 0.5 - 0.5 * cos(2 * pi * j / (B - 1))
  U <- sum(w^2)
  
  ## Fourier frequencies
  nfreq <- floor(B / 2) + 1
  freq <- (0:(nfreq - 1)) * samprate / B
  
  ## Initialize joint spectral matrix
  S <- array(0 + 0i, dim = c(nfreq, d, d))
  
  ## Welch spectral estimator
  for (s in starts) {
    
    ind <- s:(s + B - 1)
    Z.block <- Z[ind, , drop = FALSE]
    
    if (demean) {
      Z.block <- sweep(Z.block, 2, colMeans(Z.block), "-")
    }
    
    ## Apply Hann window and compute FFT
    Z.block <- Z.block * w
    D <- apply(Z.block, 2, FUN = stats::fft)
    D <- D[seq_len(nfreq), , drop = FALSE]
    
    ## Cross-periodogram matrices
    for (k in seq_len(nfreq)) {
      dk <- matrix(D[k, ], ncol = 1)
      S[k, , ] <- S[k, , ] + (dk %*% Conj(t(dk))) / (samprate * U)
    }
  }
  
  ## Average periodogram matrices across blocks
  S <- S / bnum
  
  ## Convert to one-sided spectral density
  if (B %% 2 == 0) {
    if (nfreq > 2) {
      S[2:(nfreq - 1), , ] <- 2 * S[2:(nfreq - 1), , ]
    }
  } else {
    if (nfreq > 1) {
      S[2:nfreq, , ] <- 2 * S[2:nfreq, , ]
    }
  }
  
  ## Indices for X and Y
  ix <- seq_len(p)
  iy <- p + seq_len(q)
  
  ## Inverse square root of a Hermitian matrix
  inv_sqrt_hermitian <- function(A) {
    
    A <- (A + Conj(t(A))) / 2
    
    if (ridge > 0) {
      A <- A + ridge * mean(Re(diag(A))) * diag(nrow(A))
    }
    
    ee <- eigen(A, symmetric = TRUE)
    lambda <- Re(ee$values)
    V <- ee$vectors
    
    keep <- lambda > tol * max(lambda)
    
    if (!any(keep)) {
      return(NULL)
    }
    
    Vkeep <- V[, keep, drop = FALSE]
    lkeep <- lambda[keep]
    
    Vkeep %*% diag(1 / sqrt(lkeep), nrow = length(lkeep)) %*% Conj(t(Vkeep))
  }
  
  ## Canonical coherence
  cancoh.est <- rep(NA_real_, nfreq)
  
  for (k in seq_len(nfreq)) {
    
    Sk <- S[k, , ]
    
    Sxx <- Sk[ix, ix, drop = FALSE]
    Syy <- Sk[iy, iy, drop = FALSE]
    Sxy <- Sk[ix, iy, drop = FALSE]
    
    Sxx.mhalf <- inv_sqrt_hermitian(Sxx)
    Syy.mhalf <- inv_sqrt_hermitian(Syy)
    
    if (is.null(Sxx.mhalf) || is.null(Syy.mhalf)) {
      next
    }
    
    ## Whitened cross-spectrum
    C <- Sxx.mhalf %*% Sxy %*% Syy.mhalf
    
    ## Maximum squared canonical coherence
    rho <- Re(svd(C, nu = 0, nv = 0)$d)
    rho <- pmin(pmax(rho, 0), 1)
    cancoh.est[k] <- max(rho)^2
  }
  
  return(list(freq = freq, cancoh = cancoh.est))
}



coh <- function(x, y, blength = 100) {
  
  Tlength <- length(y)
  bnum <- Tlength / blength
  nfreq <- blength / 2
  
  ## Initialize auto- and cross-periodograms
  perio_x <- matrix(NA_real_, nrow = nfreq, ncol = bnum)
  perio_y <- matrix(NA_real_, nrow = nfreq, ncol = bnum)
  perio_xy <- matrix(NA_complex_, nrow = nfreq, ncol = bnum)
  
  ## Fourier frequency indices (excluding zero)
  freq.id <- seq_len(nfreq) + 1L
  
  ## Calculate block periodograms
  for (b in seq_len(bnum)) {
    
    ind <- blength * (b - 1) + seq_len(blength)
    
    fftx <- stats::fft(x[ind])[freq.id]
    ffty <- stats::fft(y[ind])[freq.id]
    
    perio_x[, b] <- Mod(fftx)^2 / blength
    perio_y[, b] <- Mod(ffty)^2 / blength
    perio_xy[, b] <- fftx * Conj(ffty) / blength
  }
  
  ## Fourier frequencies
  freq <- seq_len(nfreq) / blength
  
  ## Average auto- and cross-periodograms
  fxx <- rowMeans(perio_x)
  fyy <- rowMeans(perio_y)
  fxy <- rowMeans(perio_xy)
  
  ## Magnitude-squared coherence
  coh.est <- Mod(fxy)^2 / (fxx * fyy)
  
  return(list(
    freq = freq,
    coh = coh.est
  ))
}




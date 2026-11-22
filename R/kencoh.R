kencoh <- function(x, y, samprate = 100,
                   fband = c("delta", "theta", "alpha", "beta", "gamma")) {
  
  Tlength <- nrow(x)
  p <- ncol(x)
  q <- ncol(y)
  
  ## Extract frequency-specific signals for X
  x.f <- array(NA_real_, dim = c(Tlength, p, length(fband)))
  
  for (ch in seq_len(p)) {
    x.band <- extract.fband(x[, ch], samprate = samprate, fband = fband)
    
    for (f in seq_along(fband)) {
      x.f[, ch, f] <- x.band[[f]]
    }
  }
  
  ## Extract frequency-specific signals for Y
  y.f <- array(NA_real_, dim = c(Tlength, q, length(fband)))
  
  for (ch in seq_len(q)) {
    y.band <- extract.fband(y[, ch], samprate = samprate, fband = fband)
    
    for (f in seq_along(fband)) {
      y.f[, ch, f] <- y.band[[f]]
    }
  }
  
  ## Estimate KenCoh for each frequency band
  kencoh.est <- rep(NA_real_, length(fband))
  
  for (f in seq_along(fband)) {
    kencoh.est[f] <- kencoh.core(
      dat = cbind(x.f[, , f], y.f[, , f]), nlags = 1, p = p, q = q
    )$kencoh
  }
  
  return(kencoh.est)
}


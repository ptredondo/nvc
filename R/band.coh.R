band.coh <- function(x, y, samprate = 100, l_max = 50,
                     fband = c("delta", "theta", "alpha", "beta", "gamma")) {
  
  Tlength <- length(x)
  
  ## Extract frequency-specific signals
  x.f <- extract_rhythms(x, samprate = samprate, fband = fband)
  y.f <- extract_rhythms(y, samprate = samprate, fband = fband)
  
  ## Define time indices and lags
  time.sub <- (l_max + 1):(Tlength - l_max)
  lags <- (-l_max):l_max
  
  ## Compute squared cross-correlations
  coh.comp <- matrix(
    NA_real_,
    nrow = length(lags),
    ncol = length(fband)
  )
  
  for (f in seq_along(fband)) {
    for (l in lags) {
      coh.comp[l + l_max + 1, f] <- stats::cor(
        x.f[[f]][time.sub],
        y.f[[f]][time.sub - l]
      )^2
    }
  }
  
  ## Maximum squared cross-correlation for each frequency band
  coh.est <- apply(coh.comp, MARGIN = 2, FUN = max, na.rm = TRUE)
  
  return(coh.est)
}



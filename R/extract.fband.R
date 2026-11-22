extract.fband <- function(x, samprate,
                            fband = c("delta", "theta", "alpha", "beta", "gamma"),
                            method = "butter", twosided = TRUE) {
  
  ## Define frequency bands in Hz
  bands <- list(
    delta  = c(0.5, 4),
    delta2 = c(0.5, 3),
    delta3 = c(1, 4),
    theta  = c(4, 8),
    theta2 = c(4, 7),
    alpha  = c(8, 12),
    alpha2 = c(9, 12),
    beta   = c(12, 30),
    beta2  = c(21, 25),
    gamma  = c(30, 45),
    gamma2 = c(36, 40)
  )
  
  ## Initialize output
  y <- vector("list", length(fband))
  
  for (i in seq_along(fband)) {
    
    ## Normalize frequencies by the Nyquist frequency
    cutoff <- bands[[fband[i]]] / (0.5 * samprate)
    
    ## Construct bandpass filter
    if (method == "butter") {
      filter_coeffs <- signal::butter(
        4, cutoff, type = "pass"
      )
    } else if (method == "fir") {
      filter_coeffs <- signal::fir1(
        100, cutoff, type = "pass"
      )
    }
    
    ## Apply filter
    if (twosided) {
      y[[i]] <- as.numeric(signal::filtfilt(filter_coeffs, x))
    } else {
      y[[i]] <- as.numeric(signal::filter(filter_coeffs, x))
    }
  }
  
  return(y)
}


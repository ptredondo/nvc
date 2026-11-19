sim.fband <- function(nsec, samprate, fband) {
  
  Tlength <- nsec * samprate
  Z <- array(NA, dim = c(Tlength, length(fband)))
  
  for (f in seq_along(fband)) {
    
    if (fband[f] == "delta") {
      
      peaklocation <- 2
      sharpness <- 0.15
      
    } else if (fband[f] == "theta") {
      
      peaklocation <- 6
      sharpness <- 0.15
      
    } else if (fband[f] == "alpha") {
      
      peaklocation <- 10
      sharpness <- 0.15
      
    } else if (fband[f] == "beta2") {
      
      peaklocation <- 23
      sharpness <- 0.15
      
    } else if (fband[f] == "gamma2") {
      
      peaklocation <- 38
      sharpness <- 0.15
      
    }
    
    psi <- peaklocation / samprate
    
    par <- c(
      2 * exp(-sharpness) * cos(2 * pi * psi),
      -exp(-2 * sharpness)
    )
    
    Z[, f] <- scale(
      extract_rhythms(
        x = stats::arima.sim(model = list(ar = par), n = Tlength),
        samprate = samprate,
        fband = fband[f]
      )[[1]]
    )
  }
  
  return(Z)
}


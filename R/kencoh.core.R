kencoh.core <- function(dat, nlags, p, q) {
  
  dat <- as.matrix(dat)
  n <- nrow(dat)
  
  idxX <- seq_len(p)
  idxY <- p + seq_len(q)
  
  ## Within-group Kendall correlation matrices
  tau.xx <- pcaPP::cor.fk(dat[, idxX, drop = FALSE])
  tau.yy <- pcaPP::cor.fk(dat[, idxY, drop = FALSE])
  
  Rxx <- sin((pi / 2) * tau.xx)
  Ryy <- sin((pi / 2) * tau.yy)
  
  ## Invert within-group matrices
  Rxx.inv <- tryCatch(solve(Rxx), error = function(e) NULL)
  Ryy.inv <- tryCatch(solve(Ryy), error = function(e) NULL)
  
  if (is.null(Rxx.inv) || is.null(Ryy.inv)) {
    return(list(kencoh = NA_real_, lag.id = NA_integer_,
                weights.y = NULL, weights.x = NULL))
  }
  
  ## Preserve original lag alignment
  rows.current <- (nlags + 1):n
  Y.fixed <- dat[rows.current, idxY, drop = FALSE]
  
  ## Initialize storage
  cc <- numeric(nlags)
  Ryx.list <- vector("list", nlags)
  
  ## Search over lags
  for (l in seq_len(nlags)) {
    
    lag <- l - 1
    X.lag <- dat[rows.current - lag, idxX, drop = FALSE]
    
    ## Kendall cross-correlation matrix
    tau.all <- pcaPP::cor.fk(cbind(X.lag, Y.fixed))
    tau.yx <- tau.all[idxY, idxX, drop = FALSE]
    
    Ryx <- sin((pi / 2) * tau.yx)
    Ryx.list[[l]] <- Ryx
    
    ## Leading eigenvalue
    val <- tryCatch(
      kencoh.eigen(Ryx, Ryy.inv, Rxx.inv, vectors = FALSE),
      error = function(e) NA_real_
    )
    
    if (is.na(val)) {
      return(list(kencoh = NA_real_, lag.id = NA_integer_,
                  weights.y = NULL, weights.x = NULL))
    }
    
    cc[l] <- Mod(val)
  }
  
  ## Find best lag
  kencoh.est <- max(cc)
  lag.id <- which(cc == kencoh.est)
  lag.best <- lag.id[1]
  
  ## Compute eigenvectors at the best lag
  best <- tryCatch(
    kencoh.eigen(Ryx.list[[lag.best]], Ryy.inv, Rxx.inv, vectors = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(best)) {
    return(list(kencoh = NA_real_, lag.id = NA_integer_,
                weights.y = NULL, weights.x = NULL))
  }
  
  return(list(
    kencoh = kencoh.est, lag.id = lag.id,
    weights.y = best$eigvec.y[, 1],
    weights.x = best$eigvec.x[, 1]
  ))
}



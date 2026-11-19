xi <- function(x, y, nboot = 0, ncores = 1) {
  n <- as.numeric(length(y))
  R <- rank(y, ties.method = "random")
  L <- n - R + 1
  M_i <- RANN::nn2(cbind(x), cbind(x), k = 2)$nn.idx[, 2]
  
  num <- n * pmin(R, R[M_i]) - L^2
  den <- L * (n - L)
  xi <- sum(num) / sum(den)
  
  if (nboot > 0) {
    
    boot.chunk <- function(B) {
      out <- numeric(B)
      
      for (b in seq_len(B)) {
        ind <- sample.int(n)
        L <- n - R[ind] + 1
        num <- n * pmin(R[ind], R[M_i]) - L^2
        den <- L * (n - L)
        out[b] <- sum(num) / sum(den)
      }
      
      return(out)
    }
    
    chunk.sizes <- rep(nboot %/% ncores, ncores)
    
    if (nboot %% ncores > 0) {
      chunk.sizes[seq_len(nboot %% ncores)] <-
        chunk.sizes[seq_len(nboot %% ncores)] + 1L
    }
    
    boot.result <- parallel::mclapply(
      chunk.sizes, boot.chunk,
      mc.cores = ncores, mc.set.seed = TRUE
    )
    
    xi_comp <- unlist(boot.result, use.names = FALSE)
    pval <- (1 + sum(xi_comp > xi)) / (nboot + 1)
    
    return(list(xi = xi, pval = pval))
  } else {
    return(xi)
  }
}



xi.vec <- function(x, y, perm.inv = TRUE, nperm = NULL,
                   nboot = 0, ncores = NULL) {
  
  ## Use double precision to avoid integer overflow in xi_M.
  n <- as.numeric(nrow(y))
  q <- ncol(y)
  
  ## Fast xi given ranks and an NN vector
  C <- n * (n + 1) * (2 * n + 1) / 6
  D <- n * (n^2 - 1) / 6
  
  xi_M <- function(R, M) {
    (n * sum(pmin(R, R[M])) - C) / D
  }
  
  ## Nearest neighbor based on general matrix z
  NN <- function(z) {
    z <- as.matrix(z)
    RANN::nn2(z, z, k = 2)$nn.idx[, 2]
  }
  
  ## Observed ranks
  R.obs <- apply(y, 2, rank, ties.method = "random")
  
  if (q == 1) {
    R.obs <- matrix(R.obs, ncol = 1)
  }
  
  ## Variable permutations
  if (perm.inv) {
    perm.mat <- gtools::permutations(n = q, r = q)
    
    if (!is.null(nperm)) {
      id <- sample(seq_len(nrow(perm.mat)),
                   min(nperm, nrow(perm.mat)),
                   replace = FALSE)
      
      perm.mat <- perm.mat[id, , drop = FALSE]
    }
  } else {
    perm.mat <- matrix(seq_len(q), nrow = 1)
  }
  
  np <- nrow(perm.mat)
  
  ## Encode each conditioning set S using a bit mask
  ## mask = 0 corresponds to S = empty set
  mask.mat <- matrix(0L, nrow = np, ncol = q)
  
  if (q >= 2) {
    for (p in seq_len(np)) {
      mask <- 0L
      
      for (i in 2:q) {
        prev <- perm.mat[p, i - 1]
        mask <- bitwOr(mask, bitwShiftL(1L, prev - 1L))
        mask.mat[p, i] <- mask
      }
    }
  }
  
  ## A state is identified by target response j + conditioning set S
  state.code <- perm.mat + q * mask.mat
  unique.states <- unique(as.vector(state.code))
  
  state.map <- matrix(match(as.vector(state.code), unique.states),
                      nrow = np, ncol = q)
  
  nstates <- length(unique.states)
  state.target <- ((unique.states - 1L) %% q) + 1L
  state.mask <- (unique.states - 1L) %/% q
  used.masks <- sort(unique(state.mask))
  
  ## Convert bit mask to variable indices
  mask_to_set <- function(mask) {
    if (mask == 0L) {
      return(integer(0))
    }
    
    which(vapply(0:(q - 1), function(k) {
      bitwAnd(mask, bitwShiftL(1L, k)) != 0L
    }, logical(1)))
  }
  
  mask.set <- lapply(used.masks, mask_to_set)
  names(mask.set) <- as.character(used.masks)
  
  ## Precompute state indices associated with each mask
  state.ids.by.mask <- split(seq_len(nstates), as.character(state.mask))
  
  ## =========================================================
  ## OBSERVED STATISTIC
  ## =========================================================
  
  num.state <- numeric(nstates)
  den.state <- numeric(nstates)
  
  ## Cache observed NN graphs
  M.num.obs <- vector("list", length(used.masks))
  M.den.obs <- vector("list", length(used.masks))
  
  names(M.num.obs) <- as.character(used.masks)
  names(M.den.obs) <- as.character(used.masks)
  
  for (mask in used.masks) {
    mask.char <- as.character(mask)
    S <- mask.set[[mask.char]]
    
    if (length(S) == 0) {
      M.num.obs[[mask.char]] <- NN(x)
    } else {
      M.num.obs[[mask.char]] <- NN(cbind(x, y[, S, drop = FALSE]))
    }
    
    if (length(S) > 0) {
      M.den.obs[[mask.char]] <- NN(y[, S, drop = FALSE])
    }
  }
  
  ## Calculate each unique observed statistic once
  for (s in seq_len(nstates)) {
    j <- state.target[s]
    mask.char <- as.character(state.mask[s])
    
    num.state[s] <- xi_M(R.obs[, j], M.num.obs[[mask.char]])
    
    if (state.mask[s] == 0L) {
      den.state[s] <- 0
    } else {
      den.state[s] <- xi_M(R.obs[, j], M.den.obs[[mask.char]])
    }
  }
  
  num.mat <- matrix(num.state[state.map], nrow = np, ncol = q)
  den.mat <- matrix(den.state[state.map], nrow = np, ncol = q)
  
  den.sum <- rowSums(den.mat)
  
  xi.all <- 1 - (q - rowSums(num.mat)) / (q - den.sum)
  xi.vec <- mean(xi.all)
  
  if (nboot <= 0) {
    return(xi.vec)
  }
  
  ## =========================================================
  ## PARALLEL PERMUTATION NULL
  ## =========================================================
  
  if (is.null(ncores)) {
    ncores <- max(1L, parallel::detectCores() - 1L)
  }
  
  ncores <- min(ncores, nboot)
  
  ## M(X) never changes
  M.X <- M.num.obs[["0"]]
  
  ## One worker handles several bootstrap replicates
  boot.chunk <- function(B) {
    
    out <- numeric(B)
    
    for (b in seq_len(B)) {
      
      ## Common row permutation preserves within-Y dependence
      ind <- sample.int(n)
      yb <- y[ind, , drop = FALSE]
      Rb <- R.obs[ind, , drop = FALSE]
      
      num.state.b <- numeric(nstates)
      
      ## One NN graph per conditioning subset
      for (mask in used.masks) {
        
        mask.char <- as.character(mask)
        S <- mask.set[[mask.char]]
        
        if (length(S) == 0) {
          M <- M.X
        } else {
          M <- NN(cbind(x, yb[, S, drop = FALSE]))
        }
        
        ## All target variables using this conditioning set
        ids <- state.ids.by.mask[[mask.char]]
        
        for (s in ids) {
          j <- state.target[s]
          num.state.b[s] <- xi_M(Rb[, j], M)
        }
      }
      
      num.mat.b <- matrix(num.state.b[state.map], nrow = np, ncol = q)
      
      xi.all.b <- 1 - (q - rowSums(num.mat.b)) / (q - den.sum)
      out[b] <- mean(xi.all.b)
    }
    
    return(out)
  }
  
  ## Split bootstrap replicates approximately equally across cores
  chunk.sizes <- rep(nboot %/% ncores, ncores)
  
  if (nboot %% ncores > 0) {
    chunk.sizes[seq_len(nboot %% ncores)] <-
      chunk.sizes[seq_len(nboot %% ncores)] + 1L
  }
  
  ## Parallel execution
  boot.result <- parallel::mclapply(
    chunk.sizes,
    boot.chunk,
    mc.cores = ncores,
    mc.set.seed = TRUE
  )
  
  xi.vec.boot <- unlist(boot.result, use.names = FALSE)
  
  pval <- (1 + sum(xi.vec.boot >= xi.vec)) / (nboot + 1)
  
  return(list(xi.vec = xi.vec, pval = pval))
}



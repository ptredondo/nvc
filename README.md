# nvc (development draft 0.1.0)

Nonlinear Vector Coherence for multivariate time series, based on the supplied
JASA ACS manuscript implementation. All 12 supplied functions are included;
10 are exported. The two KenCoh helpers are internal.

## Install locally

In R, from the directory containing the source archive:

```r
install.packages(c("RANN", "gtools", "LSTS", "signal", "pcaPP", "igraph"))
install.packages("nvc_0.1.0.tar.gz", repos = NULL, type = "source")
library(nvc)
```

If a dependency is unavailable for your R version, install a compatible
version before proceeding. To develop from the extracted `nvc/` directory,
use `devtools::load_all("nvc")` and `devtools::check("nvc")`.

## Quick start

```r
set.seed(2026)
x <- matrix(rnorm(4000), ncol = 2)
y <- matrix(rnorm(4000), ncol = 2)
fit <- nvc(x, y, blength = 100)
plot(fit$freq, fit$nvc, type = "l", xlab = "Frequency (cycles/sample)",
     ylab = "NVC")
fit.test <- nvc(x, y, blength = 100, nboot = 99, ncores = 1)
head(fit.test$pval.adj)
```

Function help is available via `?nvc`, `?xi.vec`, `?coh`, `?cancoh`,
`?kencoh`, and the other exported names.

## Development and reproducibility notes

- The uploaded code is retained, including the original calculations, defaults,
  frequency conventions, and permutation comparison rules. The package adds
  documentation, metadata, and smoke tests. No statistical algorithm has been
  silently revised.
- `nvc` and `coh` use positive normalized Fourier frequencies. `cancoh` includes
  zero and reports frequencies in Hz when `samprate` is supplied.
- For `nvc`, use matrices with at least two columns in Y and complete blocks of
  even length. `freq = NULL` uses all positive Fourier frequencies. To select
  specific bins, pass a unique, increasing subset of
  `seq(0, 0.5, length = blength / 2 + 1)[-1]` for exact matching.
  Singleton Y still loses its matrix dimension in array indexing.
- `xi` currently permutes one rank term but not its neighbor rank term in
  permutation replicates; it also uses strict `>` whereas `xi.vec` uses `>=`.
  Review this with the paper authors before relying on univariate p-values.
- `xi.vec` enumerates all response permutations before subsampling; this becomes
  costly at large q. Use small response dimensions. Row permutations require
  exchangeability under the null; assess this assumption for time-series blocks.
- Use `ncores = 1` on Windows; the existing parallel implementation uses forked
  workers on supported platforms. For reproducible tests explicitly set seeds.
- `sim.fband` is restricted to the original simulation bands delta, theta,
  alpha, beta2, gamma2. Other band names are not supported by that simulator.
- KenCoh helpers are adapted from Talento's original repository (as identified
  by the package author), with computational revisions. Confirm the precise
  repository, authors, citation, license, and permission to redistribute before
  release. Not exporting a function does not remove licensing obligations.
- Original package code is licensed under MIT (see LICENSE.md). The upstream
  KenCoh helpers are excluded from that grant pending license verification.
  Keep this development draft private until those terms and author credits
  are confirmed and the manuscript citation is finalized.
- R was not available in the packaging environment, so installation and runtime
  tests must be run locally before treating this as a validated release.

## Tests

After extracting the archive, run `R CMD build nvc` and
`R CMD check nvc_0.1.0.tar.gz` (or `devtools::check("nvc")`).
The base-R tests are in `tests/smoke.R`.

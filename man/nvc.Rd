\name{nvc}
\alias{nvc}
\alias{xi}
\alias{xi.vec}
\title{Rank-Based Dependence and Nonlinear Vector Coherence}
\description{Estimate the rank-based coefficient for a scalar response, its
multivariate extension, or the corresponding frequency-domain measure using
block periodograms. Nonlinear vector coherence (NVC), proposed by Redondo,
Huser and Ombao (2025), measures nonlinear frequency-specific dependence
between vector time series and can be used to characterize region-to-region
functional brain connectivity. Optional permutation tests assess independence.}
\usage{
xi(x, y, nboot = 0, ncores = 1)
xi.vec(x, y, perm.inv = TRUE, nperm = NULL, nboot = 0, ncores = NULL)
nvc(x, y, perm.inv = TRUE, nperm = NULL, blength = 100,
    freq = NULL, b.scale = FALSE, nboot = 0, ncores = 1)
}
\arguments{
\item{x}{Numeric predictors; rows are observations, columns are variables.}
\item{y}{Numeric response vector for \code{xi}, numeric response matrix for
\code{xi.vec}, or a multivariate time-series matrix with at least two columns
for the current \code{nvc} implementation.}
\item{perm.inv}{Average over permutations of response-variable ordering.}
\item{nperm}{Number of response orderings sampled without replacement;
\code{NULL} uses all orderings.}
\item{blength}{Positive even block length. Supply at least two complete blocks
and a series length divisible by the block length.}
\item{freq}{Frequency labels; \code{NULL} uses the complete positive Fourier
grid. Custom labels currently correspond to the first bins, not arbitrary bins.}
\item{b.scale}{Standardize each component within each block before computing periodograms.}
\item{nboot}{Number of permutation replicates; zero skips inference.}
\item{ncores}{Number of workers. For \code{xi.vec}, \code{NULL} selects available
cores minus one, subject to the replicate count. Use one on Windows.}
}
\value{\code{xi} and \code{xi.vec} return a scalar without inference; otherwise
a list containing \code{xi} or \code{xi.vec} and \code{pval}, respectively.
\code{nvc} returns a list with normalized \code{freq} and \code{nvc}; with
inference it also returns \code{pval} and BH-adjusted \code{pval.adj}.}
\details{NVC estimates frequency-specific dependence by applying the
multivariate rank-based dependence measure to componentwise block
periodograms at each Fourier frequency (Redondo et al., 2025).
When \code{perm.inv = TRUE}, the multivariate estimator averages over selected
component orderings.
The number of possible orderings grows factorially. The permutation test requires
exchangeability under the null. Negative finite-sample estimates may occur.
The univariate permutation test retains the original manuscript code and should
be reviewed before use; see the package README.}
\references{
Redondo, P. V., Huser, R., and Ombao, H. (2025).
Nonlinear Coherence for Vector Time Series: Defining Region-to-Region
Functional Brain Connectivity. \emph{arXiv preprint} arXiv:2511.14417.
\url{https://arxiv.org/abs/2511.14417}.
}
\examples{
set.seed(1)
x <- matrix(rnorm(2000), ncol = 2)
y <- matrix(rnorm(2000), ncol = 2)
xi(x, y[, 1])
xi.vec(x, y)
fit <- nvc(x, y, blength = 100)
plot(fit$freq, fit$nvc, type = "l")
}

\name{coherence}
\alias{coh}
\alias{band.coh}
\alias{cancoh}
\alias{kencoh}
\title{Comparison Measures of Spectral and Bandwise Dependence}
\description{Compute magnitude-squared coherence, maximum squared lagged
correlation of filtered signals, maximum squared canonical coherence (CanCoh),
or Kendall's tau-based canonical coherence (KenCoh).}
\usage{
coh(x, y, blength = 100)
band.coh(x, y, samprate = 100, l_max = 50,
         fband = c("delta", "theta", "alpha", "beta", "gamma"))
cancoh(x, y, samprate = 1, blength = 100, demean = TRUE,
       ridge = 0, tol = 1e-10)
kencoh(x, y, samprate = 100,
       fband = c("delta", "theta", "alpha", "beta", "gamma"))
}
\arguments{
\item{x, y}{Equal-length numeric vectors for \code{coh} and \code{band.coh};
numeric time-series matrices with equal row counts for \code{cancoh} and \code{kencoh}.}
\item{blength}{Number of observations per block; even for \code{coh}.}
\item{samprate}{Sampling rate in Hz (one gives normalized frequencies for \code{cancoh}).}
\item{l_max}{Maximum absolute lag in samples for the bandwise correlation.}
\item{fband}{Character vector of frequency bands; see \code{extract.fband}.}
\item{demean}{Remove component-specific means in each Welch block.}
\item{ridge}{Nonnegative ridge regularization factor for spectral matrices.}
\item{tol}{Relative eigenvalue threshold for spectral matrix pseudoinversion.}
}
\value{\code{coh} returns \code{freq} and \code{coh}; \code{cancoh} returns
\code{freq} and \code{cancoh}. Both bandwise functions return a numeric vector
in the order of \code{fband}.}
\details{\code{coh} uses nonoverlapping block periodograms and excludes zero.
\code{cancoh} estimates maximum squared canonical coherence (Brillinger, 1969)
using nonoverlapping Hann-tapered blocks and includes zero.
\code{band.coh} maximizes squared lagged correlation after bandpass filtering
(Ombao and Pinto, 2024). \code{kencoh} evaluates lag zero
using Kendall's tau-based canonical coherence (Talento et al., 2024).
The computations of KenCoh are implemented by adapting the original code of
Talento et al. (2024). Please contact the authors for access to the original code.}
\references{
Brillinger, D. R. (1969). The canonical analysis of stationary time series.
\emph{Multivariate Analysis}, 2, 331--350. Academic Press, New York.

Ombao, H., and Pinto, M. (2024). Spectral dependence.
\emph{Econometrics and Statistics}, 32, 122--159.

Talento, M. S. D., Roy, S., Vallejo, T. R., Ibrahim, L. A., and Ombao, H. C.
(2024). KenCoh: A ranked-based canonical coherence. \emph{arXiv preprint}
arXiv:2412.10521. \url{https://arxiv.org/abs/2412.10521}.
}
\examples{
set.seed(1)
x <- rnorm(1000)
y <- rnorm(1000)
coh(x, y)
band.coh(x, y, l_max = 2, fband = "alpha")
cancoh(cbind(x, y), cbind(y, x))
}

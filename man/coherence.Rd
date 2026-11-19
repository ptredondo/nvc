\name{coherence}
\alias{coh}
\alias{band.coh}
\alias{cancoh}
\alias{kencoh}
\title{Comparison Measures of Spectral and Bandwise Dependence}
\description{Compute magnitude-squared coherence, maximum squared lagged
correlation of filtered signals, maximum squared canonical coherence, or KenCoh.}
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
\item{fband}{Character vector of frequency bands; see \code{extract_rhythms}.}
\item{demean}{Remove component-specific means in each Welch block.}
\item{ridge}{Nonnegative ridge regularization factor for spectral matrices.}
\item{tol}{Relative eigenvalue threshold for spectral matrix pseudoinversion.}
}
\value{\code{coh} returns \code{freq} and \code{coh}; \code{cancoh} returns
\code{freq} and \code{cancoh}. Both bandwise functions return a numeric vector
in the order of \code{fband}.}
\details{\code{coh} uses nonoverlapping block periodograms and excludes zero.
\code{cancoh} uses nonoverlapping Hann-tapered blocks and includes zero.
\code{band.coh} maximizes squared lagged correlation after bandpass filtering;
it is not conventional spectral coherence. \code{kencoh} evaluates lag zero
using the revised Kendall-based implementation adapted from Talento's code.
The original citation and licensing must be confirmed before public distribution.}
\examples{
set.seed(1)
x <- rnorm(1000)
y <- rnorm(1000)
coh(x, y)
band.coh(x, y, l_max = 2, fband = "alpha")
cancoh(cbind(x, y), cbind(y, x))
}

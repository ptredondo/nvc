\name{extract.fband}
\alias{extract.fband}
\alias{sim.fband}
\title{Extract and Simulate Frequency-Band Signals}
\description{Filter a time series into named EEG frequency bands or simulate
standardized band-specific AR(2) oscillations.}
\usage{
extract.fband(x, samprate,
    fband = c("delta", "theta", "alpha", "beta", "gamma"),
    method = "butter", twosided = TRUE)
sim.fband(nsec, samprate, fband)
}
\arguments{
\item{x}{Numeric signal vector.}
\item{samprate}{Sampling rate in Hz. The Nyquist frequency must exceed the selected band upper cutoff.}
\item{fband}{Character vector of band names.}
\item{method}{Either \code{"butter"} (order 4) or \code{"fir"} (order 100).}
\item{twosided}{Use forward-backward filtering if true; otherwise one-pass filtering.}
\item{nsec}{Simulation length in seconds, with integer \code{nsec * samprate}.}
}
\details{Default bands in Hz: delta 0.5--4, theta 4--8, alpha 8--12,
beta 12--30, gamma 30--45. Alternative bands: delta2 0.5--3, delta3 1--4,
theta2 4--7, alpha2 9--12, beta2 21--25, gamma2 36--40.
The simulator supports only delta, theta, alpha, beta2 and gamma2, with
peak locations 2, 6, 10, 23 and 38 Hz, respectively.}
\value{\code{extract.fband} returns a list of filtered vectors. \code{sim.fband}
returns a matrix with one column per requested simulation band.}
\examples{
set.seed(1)
z <- sim.fband(10, 100, c("alpha", "beta2"))
extract.fband(z[, 1], 100, fband = "alpha")
}

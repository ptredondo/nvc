\name{generate.connect}
\alias{generate.connect}
\title{Draw a Circular EEG Connectivity Network}
\description{Write a circular 35-node connectivity diagram for seven brain regions
and five frequency bands to a PNG file.}
\usage{generate.connect(edges_input, plot.name, plot.title)}
\arguments{
\item{edges_input}{Data frame whose first two columns are node names (N1 through
N35), and whose other columns include \code{band}, \code{weight}, and \code{line_type}.}
\item{plot.name}{Output PNG path.}
\item{plot.title}{Title displayed on the diagram.}
}
\details{Nodes are ordered clockwise by region RF, RT, C, O, P, LT, LF, with
five nodes per region ordered delta, theta, alpha, beta, gamma. Edges are drawn
as supplied; use matching band endpoints to reproduce the manuscript diagrams.}
\value{Invisibly, the result of closing the PNG graphics device. The PNG is written to disk.}
\examples{
edges <- data.frame(from = "N1", to = "N6", band = "delta",
                    weight = 3, line_type = 1)
output <- tempfile(fileext = ".png")
generate.connect(edges, output, "Example")
unlink(output)
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/disparity_filter.R
\name{get_node_disparity}
\alias{get_node_disparity}
\title{Disparity measure of network nodes}
\usage{
get_node_disparity(net)
}
\arguments{
\item{net}{igraph; An undirected weighted network.}
}
\value{
A numeric vector of length \code{vcount(net)} with the disparity 
measure of each node.
}
\description{
Computes the disparity measure of each node in a given weighted network.
}
\examples{
# Compute the disparity of the nodes in the included US Airports network
disp <- get_node_disparity(net = air)

}
\references{
Serrano, M. A. et al. (2009) Extracting the multiscale backbone 
of complex weighted networks. \emph{PNAS} 106(16).

Barthelemy, M. et al. (2003) Spatial structure of the internet 
traffic. \emph{Physica A} 319.
}
\author{
Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
}

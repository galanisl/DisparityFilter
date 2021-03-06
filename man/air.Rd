% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/air_data.R
\docType{data}
\name{air}
\alias{air}
\title{Flights between the 500 busiest commercial airports in the US}
\format{An \code{igraph} object representing the airport network.}
\source{
\href{http://opsahl.co.uk/tnet/datasets/USairport500.txt}{Tore Opsahl's website}
}
\usage{
air
}
\description{
The airport network represents the connections between the 500 busiest commercial
airports in the United States. Two airports are linked if there was a flight
scheduled between them in 2002. The weights of a link correspond to the number of seats 
available on the scheduled flights between the airports.
}
\examples{
# Get the disparity measure of each airport
disp.air <- get_node_disparity(air)

}
\references{
Colizza V, Pastor-Satorras R, Vespignani A (2007) Reaction-diffusion processes 
and metapopulation models in heterogeneous networks. Nature Physics 3:276-282. doi:10.1038/nphys560
}
\keyword{datasets}

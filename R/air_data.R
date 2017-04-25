#'  Flights between the 500 busiest commercial airports in the US
#'
#' The airport network represents the connections between the 500 busiest commercial
#' airports in the United States. Two airports are linked if there was a flight
#' scheduled between them in 2002. The weights of a link correspond to the number of seats 
#' available on the scheduled flights between the airports.
#'
#' @docType data
#'
#' @format An \code{igraph} object representing the airport network.
#'
#' @keywords datasets
#'
#' @references Colizza V, Pastor-Satorras R, Vespignani A (2007) Reaction-diffusion processes 
#' and metapopulation models in heterogeneous networks. Nature Physics 3:276-282. doi:10.1038/nphys560
#'
#' @source \href{http://opsahl.co.uk/tnet/datasets/USairport500.txt}{Tore Opsahl's website}
#'
#' @examples
#' # Get the disparity measure of each airport
#' disp.air <- get_node_disparity(air)
#' 
"air"

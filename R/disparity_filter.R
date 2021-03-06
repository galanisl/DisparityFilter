#' Disparity measure of network nodes
#' 
#' Computes the disparity measure of each node in a given weighted network.
#' 
#' @param net igraph; An undirected weighted network.
#' 
#' @return A numeric vector of length \code{vcount(net)} with the disparity 
#' measure of each node.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Serrano, M. A. et al. (2009) Extracting the multiscale backbone 
#' of complex weighted networks. \emph{PNAS} 106(16).
#' @references Barthelemy, M. et al. (2003) Spatial structure of the internet 
#' traffic. \emph{Physica A} 319.
#' 
#' @examples
#' # Compute the disparity of the nodes in the included US Airports network
#' disp <- get_node_disparity(net = air)
#' 
#' @export
#' @import igraph
#' 
get_node_disparity <- function(net){
  N <- vcount(net)
  disp <- numeric(length = N)
  
  # Get edges incident to each node
  edg <- incident_edges(net, v = V(net))
  
  # Normalise the edge weights with respect to each node
  w <- lapply(edg, function(x) x$weight/sum(x$weight))
  
  # Compute the node disparities
  disp <- sapply(w, function(x) length(x) * sum(x^2))
  
  return(disp)
}

#' Disparity measures of the null model
#' 
#' Given a degree sequence, computes the average and standard deviation of the 
#' disparities for the null corresponding to each degree.
#' 
#' @param degrees numeric; A numeric vector containing a node degree sequence.
#' 
#' @return List with the following elements:
#' \item{mu}{Numeric vector with the null average disparities for the given 
#' degrees.}
#' \item{sd}{Numeric vector with the null standard deviation of disparities for 
#' the given degrees.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Serrano, M. A. et al. (2009) Extracting the multiscale backbone 
#' of complex weighted networks. \emph{PNAS} 106(16).
#'
#' @examples
#' # Compute the null disparities for the spectrum of node degrees of the 
#' # included US Airports network
#' null_disp <- get_node_null_disparity(degrees = seq(from = 1, 
#' to = max(degree(air)), by = 1))
#' 
#' @export
#' 
get_node_null_disparity <- function(degrees){
  mu <- 2*degrees/(degrees + 1)
  sig <- sqrt((degrees^2) * 
                ((20 + 4 * degrees)/
                   ((degrees + 1)*(degrees + 2)*(degrees + 3)) - 
                   (4/((degrees + 1)^2))
                 )
              )
  return(list(mu = mu, sd = sig))
}

#' Visualise node disparities as a function of node degree
#' 
#' Given an undirected weighted network and node disparities precomputed with 
#' function \code{get_node_disparity}, 
#' generates a ggplot comparing node disparities as a function of node degree.
#' Reference lines for null disparities, perfect node homogeneity and perfect 
#' node heterogeneity are also shown.
#' 
#' @param net igraph; An undirected weighted network.
#' @param node_disp numeric; Numeric vector of length \code{vcount(net)} with 
#' the disparity measure of each node in the network.
#' 
#' @return A ggplot comparing node disparities as a function of node degree.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Serrano, M. A. et al. (2009) Extracting the multiscale backbone 
#' of complex weighted networks. \emph{PNAS} 106(16).
#'
#' @examples
#' # Compute the disparity of the nodes in the included US Airports network
#' disp <- get_node_disparity(net = air)
#' # Plot degree vs disparities
#' p <- plot_degree_vs_disparity(net = air, node_disp = disp)
#' 
#' @export
#' @import igraph
#' @import ggplot2
#' @importFrom scales trans_breaks trans_format math_format
#' 
plot_degree_vs_disparity <- function(net, node_disp){
  null_seq <- seq(from = 1, to = max(degree(net)), by = 1)
  null_disp <- get_node_null_disparity(degrees = null_seq)
  
  df_main <- data.frame(deg = degree(net), 
                        disp = node_disp, stringsAsFactors = F)
  df_null <- data.frame(deg = null_seq, 
                        disp = null_disp$mu + 2*null_disp$sd, 
                        stringsAsFactors = F)
  
  p <- ggplot(df_main, aes_(x = ~deg, y = ~disp)) + 
    geom_point(colour = "#7fbf7b") +
    geom_hline(yintercept = 1, linetype = 2, colour = "blue") +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "blue") +
    geom_line(data = df_null, aes_(x = ~deg, y = ~disp), colour = "red") + 
    coord_cartesian(xlim = c(1, max(degree(net))), 
                    ylim = c(1, max(degree(net)))) +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format())
    ) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format())
    ) + annotation_logticks() + labs(x = "Node degree", y = "Node disparity") + 
    theme_bw() + theme(panel.grid.minor = element_blank())
  
  return(p)
}

#' Disparity filter for an undirected weighted network
#' 
#' For each weighted edge in the given network, it computes two p-values 
#' representing its local significance for the two nodes it touches.
#' This function returns only the minium p-value for the two that are computed.
#' 
#' @param net igraph; An undirected weighted network.
#' @param deg_one_pval numeric; By default, edges of degree-one nodes are 
#' assigned a p-value of 1. This can be changed through this parameter.
#' 
#' @return An igraph object with a new edge attribute:
#' \item{pval}{A numeric vector of p-values for each weighted edge of the 
#' input network.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Serrano, M. A. et al. (2009) Extracting the multiscale backbone 
#' of complex weighted networks. \emph{PNAS} 106(16).
#' 
#' @examples
#' # Get disparity p-values for the edges of the included US Airports network
#' air_with_pval <- get_edge_disparity_pvals(net = air)
#' 
#' @export
#' @import igraph
#' 
get_edge_disparity_pvals <- function(net, deg_one_pval = 1){
  
  # Get edges incident to each node
  edg <- incident_edges(net, v = V(net))
  
  # Normalise the edge weights with respect each node
  w <- lapply(edg, function(x) x$weight/sum(x$weight))
  
  # Compute p-values for each edge (from the perspective of each node), based on 
  # beta distribution with shape parameters 1 and k-1
  pvals <- lapply(w, function(x){
    if(length(x) == 1){
      return(deg_one_pval)
    }else{
      return(1 - stats::pbeta(x, shape1 = 1, shape2 = length(x) - 1))
    }})
  
  # Consider the minimum p-value from the two corresponding to each edge end
  df <- data.frame(edg = unlist(edg), pval = unlist(pvals), 
                   stringsAsFactors = F)
  pvals <- tapply(df$pval, df$edg, min)
  
  E(net)[as.numeric(names(pvals))]$pval <- pvals
  
  return(net)
}

#' Topological analysis of the disparity filter
#' 
#' After disparity p-values have been obtained, applies different p-value 
#' thresholds to the given network and computes the fraction of remaining nodes, 
#' edges and size of the largest connected component from the perspective of the 
#' original network.
#' It does the same for a filter based on weights (global filter).
#' In addition, it recommends a disparity p-value to filter the network, keeping 
#' as many nodes as possible and the most significant edges.
#' 
#' @param net igraph; The undirected weighted network to which the disparity 
#' filter was applied. It must contain a \code{pval} edge attribute.
#' @param breaks integer; The length of the sequence of values between the 
#' maximum and minimum edge p-values.
#' 
#' @return Data frame with the following columns:
#' \item{threshold}{The different p-value thresholds considered in the 
#' analysis.}
#' \item{N}{The remaining fraction of nodes at each threshold.}
#' \item{L}{The remaining fraction of links at each threshold.}
#' \item{W}{The remaining fraction of the total weight at each threshold.}
#' \item{LCC_tot}{Fraction of nodes from the original network in the largest 
#' connected component of the filtered network.}
#' \item{LCC_bb}{Fraction of nodes from the filtered network in the largest 
#' connected component of the filtered network.}
#' \item{cc}{Clustering coefficient of the filtered network at each threshold.}
#' \item{recommended}{TRUE if the threshold is the recommended one to filter 
#' the network, FALSE otherwise.}
#' \item{filter}{Whether the entry of the data frame corresponds to the 
#' Disparity or the Global filter.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Serrano, M. A. et al. (2009) Extracting the multiscale backbone 
#' of complex weighted networks. \emph{PNAS} 106(16).
#' @references Garcia-Perez, G. et al. (2016) The hidden hyperbolic geometry of 
#' international trade: World Trade Atlas 1870-2013. 
#' \emph{Scientific Reports} 6(33441).
#'
#' @examples
#' # Get disparity p-values for the edges of the included US Airports network
#' air_with_pval <- get_edge_disparity_pvals(net = air)
#' # Analyse the topology of the networks resulting from the application of 
#' # different disparity filters
#' analysis <- analyse_disparity_filter(net = air_with_pval, breaks = 100)
#' 
#' @export
#' @import igraph
#' 
analyse_disparity_filter <- function(net, breaks = 100){
  
  disparity_cuts <- seq(from = 1, to = 0, length.out = breaks)
  
  #The weight breaks are logarithmically spaced
  global_cuts <- exp(log(10) * seq(log10(min(E(net)$weight)), 
                                   log10(max(E(net)$weight)), 
                                   length.out = breaks))
  
  remaining <- data.frame(threshold = c(disparity_cuts, global_cuts), 
                          N = numeric(length = 2*length(disparity_cuts)), 
                          L = numeric(length = 2*length(disparity_cuts)), 
                          W = numeric(length = 2*length(disparity_cuts)), 
                          LCC_tot = numeric(length = 2*length(disparity_cuts)),
                          LCC_bb = numeric(length = 2*length(disparity_cuts)), 
                          cc = numeric(length = 2*length(disparity_cuts)),
                          recommended = logical(length = 
                                                  2*length(disparity_cuts)), 
                          filter = factor(rep(c("Disparity", "Global"), 
                                              each = length(disparity_cuts)),
                                          levels = c("Disparity", "Global"), 
                                          ordered = T))
  
  for(i in 1:nrow(remaining)){
    if(remaining$filter[i] == "Disparity"){
      g <- delete_edges(net, 
                        edges = E(net)[E(net)$pval > remaining$threshold[i]])
    }else{
      g <- delete_edges(net, 
                        edges = E(net)[E(net)$weight < remaining$threshold[i]])
    }
    g <- delete_vertices(g, v = which(degree(g) < 1))
    
    remaining$N[i] <- vcount(g)/vcount(net)
    remaining$L[i] <- ecount(g)/ecount(net)
    remaining$W[i] <- sum(E(g)$weight)/sum(E(net)$weight)
    remaining$LCC_tot[i] <- ifelse(length(clusters(g)$csize) > 0, 
                                   max(clusters(g)$csize)/vcount(net), 0)
    remaining$LCC_bb[i] <- ifelse(length(clusters(g)$csize) > 0, 
                                  max(clusters(g)$csize)/vcount(g), 0)
    remaining$cc[i] <- ifelse(is.na(transitivity(g, type = "localaverage")), 
                              0, transitivity(g, type = "localaverage"))
  }
  
  # Compute the vertical distance between every point in the Lbb/Ltot<->Nbb/Ntot 
  # curve and the diagonal (where Lbb/Ltot == Nbb/Ntot)
  diag_vals <- remaining$L[remaining$filter == "Disparity"]
  d <- sqrt((remaining$L[remaining$filter == "Disparity"] - diag_vals)^2 + 
              (remaining$N[remaining$filter == "Disparity"] - diag_vals)^2)
  # The Lbb/Ltot<->Nbb/Ntot pair that maximises the vertical distance to the 
  # diagonal is set as the recommended threshold
  remaining$recommended[which.max(d)] <- TRUE
  
  return(remaining)
}

#' Visual analysis of the disparity filter applied to an undirected weighted 
#' network
#' 
#' After the application of function \code{analyse_disparity_filter} to an 
#' undirected weighted network, it generates three ggplots showing how the 
#' remaining fraction of nodes changes as a function of the remaining fraction 
#' of links and total weight as more stringent filters are applied to the 
#' network.
#' The third plot shows how the fraction of nodes in the largest connected 
#' component of the network changes as a function of the different filters 
#' applied.
#' All three plots indicate the final situation of the network for the 
#' recommended disparity filter.
#' 
#' @param disp_analysis data frame; The data frame resulting from the 
#' application of function \code{analyse_disparity_filter}.
#' 
#' @return A list of five ggplots:
#' \item{LvsN}{The remaining fraction of nodes as a function of the remaining 
#' fraction of links.}
#' \item{WvsN}{The remaining fraction of nodes as a function of the remaining 
#' fraction of total weight.}
#' \item{AvsLCC_tot}{Fraction of nodes from the original network in the LCC of 
#' the filtered network as a function of different p-value thresholds.}
#' \item{AvsLCC_bb}{Fraction of nodes from the filtered network in the LCC of 
#' the filtered network as a function of different p-value thresholds.}
#' \item{AvsCC}{Clustering coefficient as a function of different p-value 
#' thresholds.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @references Serrano, M. A. et al. (2009) Extracting the multiscale backbone 
#' of complex weighted networks. \emph{PNAS} 106(16).
#' @references Garcia-Perez, G. et al. (2016) The hidden hyperbolic geometry of 
#' international trade: World Trade Atlas 1870-2013. 
#' \emph{Scientific Reports} 6(33441).
#'
#' @examples
#' # Get disparity p-values for the edges of the included US Airports network
#' air_with_pval <- get_edge_disparity_pvals(net = air)
#' # Analyse the topology of the networks resulting from the application of 
#' # different disparity filters
#' analysis <- analyse_disparity_filter(net = air_with_pval, breaks = 100)
#' # Plot the results of the analysis
#' p <- plot_disparity_filter_analysis(disp_analysis = analysis)
#' 
#' @export
#' @import ggplot2
#' 
plot_disparity_filter_analysis <- function(disp_analysis){
  p <- vector(mode = "list", length = 5)
  names(p) <- c("LvsN", "WvsN", "AvsLCC_tot", "AvsLCC_bb", "AvsCC")
  
  p$LvsN <- ggplot(disp_analysis, 
                   aes_(~L, ~N, colour = ~filter, shape = ~filter)) + 
    geom_point(size = 4) + geom_line() +
    scale_colour_manual(values = c("#7fbf7b", "#af8dc3")) +
    geom_abline(slope = -1, intercept = 0, linetype = 2, colour = "blue") +
    geom_vline(xintercept = disp_analysis$L[disp_analysis$recommended], 
               linetype = 2, colour = "red") + 
    geom_hline(yintercept = disp_analysis$N[disp_analysis$recommended], 
               linetype = 2, colour = "red") + scale_x_reverse() + 
    labs(x = expression(L[bb]/L[tot]), y = expression(N[bb]/N[tot])) + 
    theme_bw() + 
    theme(legend.title = element_blank(), legend.background = element_blank(), 
          legend.justification = c(0,0), legend.position = c(0,0))
  
  p$WvsN <- ggplot(disp_analysis, 
                   aes_(~W, ~N, colour = ~filter, shape = ~filter)) + 
    geom_point(size = 4) + geom_line() +
    scale_colour_manual(values = c("#7fbf7b", "#af8dc3")) +
    geom_vline(xintercept = disp_analysis$W[disp_analysis$recommended], 
               linetype = 2, colour = "red") + 
    geom_hline(yintercept = disp_analysis$N[disp_analysis$recommended], 
               linetype = 2, colour = "red") + scale_x_reverse() + 
    labs(x = expression(W[bb]/W[tot]), y = expression(N[bb]/N[tot])) + 
    theme_bw() + 
    theme(legend.title = element_blank(), legend.background = element_blank(), 
          legend.justification = c(0,0), legend.position = c(0,0))
  
  #Normalise the weight for the Global filter
  disp_analysis$threshold[disp_analysis$filter == "Global"] <- 
    (disp_analysis$threshold[disp_analysis$filter == "Global"] - 
       min(disp_analysis$threshold[disp_analysis$filter == "Global"])) /
    (max(disp_analysis$threshold[disp_analysis$filter == "Global"]) - 
       min(disp_analysis$threshold[disp_analysis$filter == "Global"]))
  
  p$AvsLCC_tot <- ggplot(disp_analysis, 
                         aes_(~threshold, ~LCC_tot, colour = ~filter, 
                              shape = ~filter)) + 
    geom_point(size = 4) + geom_line() + 
    scale_colour_manual(values = c("#7fbf7b", "#af8dc3")) +
    geom_vline(xintercept = disp_analysis$threshold[disp_analysis$recommended], 
               linetype = 2, colour = "red") + 
    geom_hline(yintercept = disp_analysis$LCC_tot[disp_analysis$recommended], 
               linetype = 2, colour = "red") + scale_x_reverse() + 
    labs(x = "Disparity p-values/Normalised weight", 
         y = expression(paste("Fraction of ", N[tot], " in ", LCC[bb]))) + 
    theme_bw() + 
    theme(legend.title = element_blank(), legend.background = element_blank(), 
          legend.justification = c(0,0), legend.position = c(0,0.5))
  
  p$AvsLCC_bb <- ggplot(disp_analysis, 
                        aes_(~threshold, ~LCC_bb, colour = ~filter, 
                             shape = ~filter)) + 
    geom_point(size = 4) + geom_line() +
    scale_colour_manual(values = c("#7fbf7b", "#af8dc3")) +
    geom_vline(xintercept = disp_analysis$threshold[disp_analysis$recommended],
               linetype = 2, colour = "red") + 
    geom_hline(yintercept = disp_analysis$LCC_bb[disp_analysis$recommended], 
               linetype = 2, colour = "red") + scale_x_reverse() + 
    labs(x = "Disparity p-values/Normalised weight", 
         y = expression(paste("Fraction of ", N[bb], " in ", LCC[bb]))) + 
    theme_bw() + 
    theme(legend.title = element_blank(), legend.background = element_blank(), 
          legend.justification = c(0,0), legend.position = c(0,0.5))
  
  p$AvsCC <- ggplot(disp_analysis, 
                    aes_(~threshold, ~cc, colour = ~filter, shape = ~filter)) + 
    geom_point(size = 4) + geom_line() +
    scale_colour_manual(values = c("#7fbf7b", "#af8dc3")) +
    geom_vline(xintercept = disp_analysis$threshold[disp_analysis$recommended], 
               linetype = 2, colour = "red") + 
    geom_hline(yintercept = disp_analysis$cc[disp_analysis$recommended], 
               linetype = 2, colour = "red") + scale_x_reverse() + 
    labs(x = "Disparity p-values/Normalised weight", 
         y = "Clustering coefficient") + theme_bw() + 
    theme(legend.title = element_blank(), legend.background = element_blank(), 
          legend.justification = c(0,0), legend.position = c(0,0.5))
  
  return(p)
}

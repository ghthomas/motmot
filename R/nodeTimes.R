#' @title Get times for nodes and tips
#' @description Produces branching and tip times for ultrametric and non-ultrametric trees
#' @note \code{nodeTimes} is essentially a re-written version of the "ape" \code{branching.times}.
#' @param phy An object of class "phylo" (see ape package).
#' @return Returns a matrix corresponging the phy "edge" matrix showning internal and external node times
#' @author Mark Puttick, Emmanuel Paradis
#' @export

nodeTimes <- function(phy) {
	
	phy <- reorder(phy)
	depBranches <- node.depth.edgelength(phy)
	all <- cbind(depBranches[phy$edge[,1]], depBranches[phy$edge[,2]])
	return(max(depBranches) - all)
	
}
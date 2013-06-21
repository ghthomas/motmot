node.descendents <- function (x, phy, tip.labels = FALSE) 
{
    
	if (x <= Ntip(phy)) { return(numeric()) } else {
		node <- x
		Ntip <- length(phy$tip.label)
		ROOT <- Ntip + 1
		Nedge <- dim(phy$edge)[1]
		wbl <- !is.null(phy$edge.length)
		
		phy <- reorder(phy)
		root.node <- which(phy$edge[, 2] == node)
		start <- root.node + 1
		anc <- phy$edge[root.node, 1]
		next.anc <- which(phy$edge[-(1:start), 1] <= anc)
		keep <- if (length(next.anc)) 
        start + 0:(next.anc[1] - 1)
		else start:Nedge
		
		nodes <- phy$edge[keep, 2]    
		nodes <- unique(nodes)
		if (tip.labels == TRUE) {
			nodesTips <- vector(mode = "list", length = 2)
			nodesTips[[1]] <- nodes[nodes > Ntip]
			nodesTips[[2]] <- with(phy, tip.label[nodes[nodes <= 
								   Ntip]])
			return(nodesTips)
		}
		else {
			return(nodes)
		}}
}
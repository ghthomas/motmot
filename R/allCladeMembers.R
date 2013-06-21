allCladeMembers <- function(phy) {
	
	nodeIDs <- c(1:Ntip(phy), (Ntip(phy)+2):(length(phy$edge.length)+1)) ######## Problem here
	k <- length(nodeIDs)
    cladeMembersMatrix <- matrix(NA, ncol = k, nrow = length(phy$edge[, 
															 1]))
    for (i in 1:k) {
        nodeShiftID <- c(nodeIDs[i], node.descendents(x = nodeIDs[i], 
													  phy = phy))
        cladeMembersMatrix[, i] <- as.numeric(phy$edge[, 2] %in% nodeShiftID)
    }
	return(cladeMembersMatrix)
}
cladeIdentity <- function (phy, nodeIDs, cladeMembersObj=NULL) 
{
	
	k <- length(nodeIDs)
	
    if(is.null(cladeMembersObj)) {
		cladeMembers <- matrix(NA, ncol = k, nrow = length(phy$edge[, 
														   1]))
		for (i in 1:k) {
			nodeShiftID <- c(nodeIDs[i], node.descendents(x = nodeIDs[i], phy = phy))
			cladeMembers[, i] <- as.numeric(phy$edge[, 2] %in% nodeShiftID)
		}
	}
	
	if (is.null(cladeMembersObj)==FALSE) {
		allNodes <- c(1:Ntip(phy), (Ntip(phy)+2):(length(phy$edge.length)+1))   #######
		cladeMembers <- as.matrix(cladeMembersObj[,match(nodeIDs, allNodes)])
	}
	
    originalOrder <- colSums(cladeMembers)
    richnessOrder <- sort(originalOrder, decreasing = FALSE, 
						  index.return = TRUE)
    cladeMembersOrdered <- matrix(cladeMembers[, richnessOrder$ix], 
								  ncol = length(nodeIDs))
	
	if (k>1) {
		for (i in 2:k){
			if (i ==2 ) {	cladeMembersOrdered[,i] <- cladeMembersOrdered[,i] - cladeMembersOrdered[,1:(i-1)]}
			else {cladeMembersOrdered[,i] <- cladeMembersOrdered[,i] - rowSums(cladeMembersOrdered[,1:(i-1)])}
		}
	}
	
    cladeMembers <- cladeMembersOrdered[, sort(richnessOrder$ix, index.return = TRUE)$ix]
    cladeMembers <- matrix(cladeMembers, ncol = length(nodeIDs))
    return(cladeMembers)
}


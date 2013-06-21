fairProportions <- function (phy, nodeCount=FALSE) {
	
	treeMatrix <- clade.matrix(phy)
	
	fpEdgeVals <- treeMatrix$edge.length / apply(treeMatrix$clade.matrix, 1, sum) 
	
# fpEdgeMatrix <- fpEdgeVals * treeMatrix$clade.matrix
	
	fpTips <- as.matrix(apply(fpEdgeVals * treeMatrix$clade.matrix, 2, sum))
	
	
	
	if (nodeCount==TRUE) {	nodeCount <- apply(treeMatrix$clade.matrix, 2, sum) - 1
		fpTips <- cbind(fpTips, nodeCount)
		colnames(fpTips) <- c("FP", "NodeCount")
		rownames(fpTips) <- phy$tip.label }
	else { rownames(fpTips) <- phy$tip.label }
	
	
	rm(treeMatrix)
	return(fpTips)
	
}
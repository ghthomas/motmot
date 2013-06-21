pic.motmot <- function (x, phy) {
out <- pic(x, phy, rescaled.tree=TRUE, var.contrasts=TRUE, scaled=FALSE)
   contr <- out[[1]]	
	nb.tip <- length(phy$tip.label)
	idx <- which(out[[2]]$edge[,1] == (nb.tip+1) )
	
	root.v <- out[[2]]$edge.length[idx]
	V <- prod(root.v)/(sum(root.v))
  
  return(list(contr = contr, root.v = root.v, V = V))
}

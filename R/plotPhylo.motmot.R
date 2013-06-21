plotPhylo.motmot <- function(phy, x=NULL, traitMedusaObject=NULL, reconType="rates", type = "phylogram", use.edge.length = TRUE, 
    show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black", 
    edge.width = 1, edge.lty = 1, font = 3, cex = par("cex"), 
    adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE, 
    label.offset = 0.5, underscore = FALSE, x.lim = NULL, y.lim = NULL, 
    direction = "rightwards", lab4ut = "horizontal", tip.color = "black", palette="hotspot.colors")
    
    
{
   	
	cool.colors <- function (n, alpha = 1) {rainbow(n, start = 3/6, end = 4/6, alpha = alpha)}
	heat.colors <- function (n, alpha = 1) { rainbow(n, start = 0, end = 1/6, alpha = alpha)}
	hotspot.colors <- function (n, alpha = 1) { rainbow(n, start = 0, end = 4/6, alpha = alpha)}  
 	
switch (reconType,
		
	"picReconstruction" = {
		if (is.vector(x)==FALSE) {stop("trait must be a vector with taxon names") }
   	
		a <- ace(x, phy, method="pic")
		reconTrait <- rep(NA, length(phy$edge[,2]))

		# internal nodes   
		reconTrait[match(names(a$ace)[2:length(a$ace)], phy$edge[,2])] <- a$ace[2:length(a$ace)]

		# tips
		idx <- match(names(x), phy$tip.label)
		reconTrait[match(idx, phy$edge[,2])] <- x

		ncolours <- Ntip(phy) + Nnode(phy) -1
		edgeColours <- rep(NA, ncolours)
    
		if (palette=="hotspot.colors"){    use.colors <- rev(hotspot.colors(ncolours))}
		if (palette=="heat.colors"){    use.colors <- rev(heat.colors(ncolours))}
		if (palette=="cool.colors"){    use.colors <- rev(cool.colors(ncolours))}
		if (palette=="combi.colors"){    use.colors <- c(cool.colors(ncolours/2), rev(heat.colors(ncolours/2)))}
		
 
		binsize <- (max(reconTrait) - min(reconTrait)) / length(reconTrait)
		color.bin <- matrix(NA, ncol=length(reconTrait), nrow=2)
		rownames(color.bin) <- c("min", "max")
 
		for (i in 1:length(reconTrait)) { 	
			color.bin[1,i] <- min(reconTrait) + (i-1) * binsize
			color.bin[2,i] <- min(reconTrait) + i* binsize
			}
 
		color.bin[2,length(reconTrait)] <- max(reconTrait)	# necessary due to rounding problems
 	
 	
		for (i in 1: ncolours) {
			color.trait <- (reconTrait[i] >= color.bin["min",]) & (reconTrait[i] <= color.bin["max",])  
			edgeColours[i] <- use.colors[which(color.trait==TRUE)]
			}
		},
		
	"rates" = {
		
		cladeRates <- traitMedusaObject$Rates[,3]
		nodes <- traitMedusaObject$Rates[,1]
		rateType <- traitMedusaObject$Rates[,2]
		
		ncolours <- Ntip(phy) + Nnode(phy) -1
		
		nodeBranch <- data.frame(as.numeric(nodes), rateType)
		
		nodes <- nodeBranch[nodeBranch[,"rateType"]=="clade",1]
		branches <- nodeBranch[nodeBranch[,"rateType"]=="branch",1]
		
   	   	cladeMembers <- cladeIdentity(phy, nodes)
		edgeColours <- rep("black", length(phy$edge[,1]))
		
		lgcladeRates <- log(as.numeric(cladeRates))
		lgcladeRatesNodes <- lgcladeRates[nodeBranch[,"rateType"]=="clade"]
		lgcladeRatesBranches <- lgcladeRates[nodeBranch[,"rateType"]=="branch"]
		
		if (palette=="hotspot.colors"){    use.colors <- rev(hotspot.colors(ncolours))}
		if (palette=="heat.colors"){    use.colors <- rev(heat.colors(ncolours))}
		if (palette=="cool.colors"){    use.colors <- rev(cool.colors(ncolours))}
		if (palette=="combi.colors"){    use.colors <- c(cool.colors(ncolours/2), rev(heat.colors(ncolours/2)))}
		
		binsize <- (max(abs(lgcladeRates)) - -max(abs(lgcladeRates))) / length(edgeColours)
		color.bin <- matrix(NA, ncol=length(edgeColours), nrow=2)
		
		rownames(color.bin) <- c("min", "max")
		
		for (i in 1:length(edgeColours)) { 	
			color.bin[1,i] <- -max(abs(lgcladeRates)) + (i-1) * binsize
			color.bin[2,i] <- -max(abs(lgcladeRates)) + i* binsize
			}
		
		color.bin[2,length(edgeColours)] <- max(abs(lgcladeRates))	# necessary due to rounding problems
		color.bin[1,1] <- -max(abs(lgcladeRates))	# necessary due to rounding problems
		
		for (i in 1:length(nodes)) {
			color.trait <- (lgcladeRatesNodes[i] >= color.bin["min",]) & (lgcladeRatesNodes[i] <= color.bin["max",]) 
			edgeColours[cladeMembers[,i]==1] <- use.colors[which(color.trait==TRUE)]
			}
		
		for (i in 1:length(branches)) {
			color.trait <- (lgcladeRatesBranches[i] >= color.bin["min",]) & (lgcladeRatesBranches[i] <= color.bin["max",]) 
			edgeColours[which(phy$edge[,2]==branches[i])] <- use.colors[which(color.trait==TRUE)]
			}
		
		}
		
		)
		
		
	plot.phylo(x=phy, edge.color = edgeColours, edge.width = edge.width, type=type, use.edge.length=use.edge.length, cex=cex, label.offset=label.offset, show.tip.label=show.tip.label, no.margin=no.margin, direction=direction, adj=adj)
	return(as.data.frame(t(rbind(use.colors, exp(color.bin)))))
    }
    
   
   
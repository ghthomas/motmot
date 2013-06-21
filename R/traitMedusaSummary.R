traitMedusaSummary <- function (traitMedusaObject=NULL, cutoff=4, AICc=TRUE, lowerBound=1e-7, upperBound=200) {
		
	y <- traitMedusaObject[[3]]
	phy <- traitMedusaObject[[4]]
	breaks <- numeric()

    if (AICc==TRUE) { datCol = 6 } else { datCol = 5 }
    
	bestModel <- traitMedusaObject[[1]][1,]
	rateType <- traitMedusaObject[[1]][2:nrow(traitMedusaObject[[1]]),2]
	
	for (i in 2:dim(traitMedusaObject[[1]])[1]) {
        if ((traitMedusaObject[[1]][i-1, datCol] - traitMedusaObject[[1]][i, datCol]) < cutoff) 
            break
			bestModel <- traitMedusaObject[[1]][i,]
			}
        
    bestModelOut <- bestModel[,which(is.na(bestModel)==FALSE)]
	bestModelOut <- bestModelOut[,2:ncol(bestModelOut)]
	
	
	
	out <- vector(mode="list", length=3)
	names(out) <- c("ModelFit", "Rates", "optimalTree")
	
	
	
	if (bestModel$node==0) { out$optimalTree <- phy 
							out$ModelFit <- bestModelOut[,2:5]
							out$Rates <- "Single rate"
							}

	
	else {
	
	
		foo <- function(param) {
			ll <- transformPhylo.ll(y, phyClade, model="clade", nodeIDs=SingleNode, cladeRates=param, rateType=whichRateType)$logLikelihood
			return(ll - bestModelOut$ML + 1.92)
		}
		
	
		nodeIDs <- as.numeric(names(bestModelOut)[6:ncol(bestModelOut)])
		cladeRates <- as.numeric(bestModelOut[,6:ncol(bestModelOut)])
		rateType <- traitMedusaObject[[1]][2:(length(cladeRates)+1),2]
		optimalTree <- transformPhylo(phy, model="clade", nodeIDs=nodeIDs, cladeRates=cladeRates, rateType=rateType)

		out$Rates <- matrix(NA, length(nodeIDs), 5, byrow=TRUE)
		out$Rates <- as.data.frame(out$Rates)
		colnames(out$Rates) <- c("node", "shiftPos", "MLRate", "LowerCI", "UpperCI")
		
		out$ModelFit <- bestModelOut[,2:5]
	
		out$optimalTree <- optimalTree

		
		for (i in 1:length(nodeIDs)) {
			
			SingleNode <- nodeIDs[i]
			whichRateType <- rateType[i]
			phyClade <- transformPhylo(optimalTree, model="clade", nodeIDs=SingleNode, cladeRates=1/cladeRates[i], rateType=whichRateType)
			LCI <- NULL
			UCI <- NULL
			
			if(foo(lowerBound) < 0) { 
				LCI <- uniroot(foo, interval = c(lowerBound, cladeRates[[i]]))$root 
			} else { LCI <- NA }
			if(foo(upperBound) < 0) {
				UCI <- uniroot(foo, interval = c(cladeRates[[i]], upperBound))$root
			} else { UCI <- NA }
			

			
			out$Rates[i,] <- aaa <- c(nodeIDs[i], rateType[i], cladeRates[i], LCI, UCI)
		}
		
		if (any(is.na(out$Rates[,3:4]))) { warning("Confidence limits fall outside the current parameter bounds - consider changing lowerBound and/or upperBound")}
												   

		
		
		 }
	return(out)
}
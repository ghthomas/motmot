transformRateMatrix <- function(rateData, rate=NULL) {
			V <- rateData$Vmat
			
			if(is.null(rate))  { rate <- c(rep(1,length(rateData$Vmat))) } else { rate <- rate }	
			
			nV <- length(rate)
			
			if(length(rate) != length(rateData$Vmat)){stop("The number of rates defined differs from the number of rate matrices")}
			
			v1 <- V[[1]]		
			rateMats <- vector(mode="list", length = nV)
			retMat <- matrix(0, nrow = dim(v1)[1], ncol = dim(v1)[2])
			
			for(i in 1:nV) {
			   rateMats[[i]] <- rate[i] * V[[i]]  
			   retMat <- retMat + rateMats[[i]]
					}
	
	return(retMat)}

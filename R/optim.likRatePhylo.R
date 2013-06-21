optim.likRatePhylo <-
function(rateData, rate=NULL, fixed = NULL, rateMIN = 0.001, rateMAX = 50, common.mean=FALSE, lambda.est=TRUE, meserr=FALSE) {
		
		if(is.null(rate))  { rate <- c(rep(1,length(rateData$Vmat))) } else { rate <- rate }	
	
		if(is.null(fixed))  { op <- fixed <- c(TRUE, rep(FALSE,length(rateData$Vmat) - 1)) } else { op <- fixed }
				
		lower <- c(rep(rateMIN, length(rateData$Vmat)), 0.00001)
		upper <- c(rep(rateMAX, length(rateData$Vmat)), 1)
	
		mvl <- make.likRatePhylo(rateData, fixed=op, common.mean=common.mean, lambda.est=lambda.est, meserr=meserr)

			vo <- try(optim(c(rate, 1, 1), mvl, method = "L-BFGS-B", lower = lower, upper = upper))
			MLRate <- vo$par[1:length(fixed)]
			Lambda <- vo$par[1+length(fixed)]
			
			fixed[which(fixed==FALSE)] <- MLRate[which(fixed==FALSE)]
			MLRate <- fixed

			ML <- -vo$value
			convergence <- vo$convergence
			n <- length(rateData$y)
			
				if(length(op)!=length(which(op==FALSE))) {
					if(common.mean==TRUE) {k <- 2 + length(which(op==FALSE) + lambda.est + meserr)
							} else { k <- (length(which(op==FALSE)) +1 + length(op)) + lambda.est + meserr }
							
							} else {
					
					if(common.mean==TRUE) {k <- 1 + length(which(op==FALSE) + lambda.est + meserr)
							} else { k <- (length(which(op==FALSE)) + length(op)) + lambda.est + meserr}
							}
							
			aic <- -2 * ML + 2 * k
			aicc <- -2 * ML + 2 * k + ((2*k*(k+1))/(n-k-1))
			
			
	ML.RatePhylo <- list(MLRate = MLRate, Lambda = Lambda, Max.lik = ML, aic = aic, aicc = aicc, convergence=convergence, n.parameters = k)
	return(ML.RatePhylo)
	
}


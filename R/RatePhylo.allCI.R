RatePhylo.allCI <-
function(rateData, MLrate=NULL, fixed=NULL, rateMIN = 0.001, rateMAX = 50, common.mean=FALSE, lambda.est=TRUE) {

	if(is.null(MLrate))  { MLrate <- c(rep(1,length(rateData$Vmat))) } else { MLrate <- MLrate }	
	if(is.null(fixed))  { fixed <- c(rep(FALSE,length(rateData$Vmat) - 1), FALSE) } else { op <- fixed }	

	
	
	n.rate <- length(MLrate)	
	all.CIs <- matrix(nrow=n.rate, ncol=2, dimnames=list(c(1:n.rate), c("Lci", "Uci")))
	
	for(i in 1:n.rate) {
			fix.now <- c(rep(TRUE, n.rate))
			fix.now[i] <- FALSE
	
		if(fixed[i] == FALSE) { 
				CI <- RatePhylo.CI(rateData, MLrate, fix.now, common.mean=common.mean, lambda.est=lambda.est)
				all.CIs[i,1] <- CI[1]
				all.CIs[i,2] <- CI[2]
				} else { all.CIs[i,1] <- NA
					all.CIs[i,2] <- NA }
					}
	
	return(all.CIs)
	}


RatePhylo.CI <-
function(rateData, MLrate=NULL, fixed, rateMIN = 0.001, rateMAX = 50, common.mean=FALSE, lambda.est=TRUE) {

	if(is.null(MLrate))  { MLrate <- c(rep(1,length(rateData$Vmat))) } else { MLrate <- MLrate }	
	
	MLrate <- as.numeric(format(MLrate))

	ML <- likRatePhylo(rateData, MLrate, common.mean=common.mean, lambda.est=lambda.est)$ll
	
		fixed.rates <- data.frame(MLrate, fixed)
		use.rate.ind <- which(fixed.rates$fixed==FALSE)
		use.rate <- fixed.rates$MLrate[use.rate.ind]
	
		var.fun <- function(vary.rate) {
		
				fixed.rates$MLrate[use.rate.ind] <- vary.rate
				test.rate <- fixed.rates$MLrate
				ll <- likRatePhylo(rateData, test.rate, common.mean=common.mean, lambda.est=lambda.est)$ll
		
		return( ll - ML + 1.92)
		 }
		 
	if(var.fun(rateMIN) < 0) { 
			rateLci <- uniroot(var.fun, interval = c(rateMIN, use.rate))$root 
			}
	if(var.fun(rateMAX) < 0) {
			rateUci <- uniroot(var.fun, interval = c(use.rate, rateMAX))$root
			}
			return(c(rateLci=rateLci, rateUci=rateUci))
			}


ML.RatePhylo <-
function(rateData, rate=NULL, fixed = NULL, pretty = TRUE, rateMIN = 0.001, rateMAX = 50, common.mean=FALSE, lambda.est=TRUE, est.CI=FALSE, meserr=FALSE, file=NULL) {
						
		if(is.null(rate))  { rate <- c(rep(1,length(rateData$Vmat))) } else { rate <- rate }	
	
		if(is.null(fixed))  { op <- fixed <- c(rep(FALSE,length(rateData$Vmat) - 1), TRUE) } else { op <- fixed }		
		
				ovp <- optim.likRatePhylo(rateData, rate, fixed, rateMIN, rateMAX, common.mean=common.mean, lambda.est=lambda.est, meserr=meserr)
				max.lik.rate <- ovp$MLRate
				lambda <- ovp$Lambda
				max.lik <- ovp$Max.lik
				singlerate <- optim.likRatePhylo(rateData, rate, rep(TRUE, length(rateData$Vmat)), rateMIN, rateMAX, common.mean=common.mean, lambda.est=lambda.est, meserr=meserr)
				lik1 <- singlerate$Max.lik
				lambda_single <- singlerate$Lambda
				D <- 2 * (max.lik - lik1)
	
	
				mlparams <- likRatePhylo(rateData, rate=max.lik.rate, common.mean=common.mean, lambda.est=lambda)
	
	
				if (common.mean==TRUE) { mu <- mlparams$mu } else {
					
					mu <- rep(NA, length(mlparams$mu))
					mu[1] <- mlparams$mu[1]
					for (i in 2:length(mlparams$mu)) { mu[i] <- mlparams$mu[1] + mlparams$mu[i] }
				}
		
														
		
				
				if (length(op) == length(which(op==FALSE)))  { k2 <- length(op) - 1 
						} else { k2 <- length(which(op==FALSE)) }
				
				if(length(op)!=length(which(op==FALSE))) {
					if(common.mean==TRUE) {k <- 2 + length(which(op==FALSE)) + lambda.est + meserr
							} else { k <- (length(which(op==FALSE)) +1 + length(op)) + lambda.est + meserr}
							
							} else {
					
					if(common.mean==TRUE) {k <- 1 + length(which(op==FALSE)) + lambda.est + meserr
							} else { k <- (length(which(op==FALSE)) + length(op)) + lambda.est + meserr}
							}
								
				pval <- 1- pchisq(D, k2)
				n <- length(rateData$y)
				
				
	if (est.CI) { CIs.rate <- RatePhylo.allCI(rateData, max.lik.rate, fixed=rep("FALSE", length(rateData$Vmat)), common.mean=common.mean)
	} else { CIs.rate <- matrix(NA, ncol=2, nrow=1) }
				aic <- -2 * max.lik + 2 * k
				aicc <- -2 * max.lik + 2 * k + ((2*k*(k+1))/(n-k-1))
				
				if(common.mean==TRUE) {
						aic.rate1 <- -2 * lik1 + 2 * 2
						aicc.rate1 <- -2 * lik1 + 2 * 2 + ((2*2*(2+1))/(n-2-1))
					} else {
						aic.rate1 <- -2 * lik1 + 2 * (1 + length(op))
						aicc.rate1 <- -2 * lik1 + 2 * (1 + length(op)) + ((2*(1 + length(op))*((1 + length(op))+1))/(n-(1 + length(op))-1))
						}
				
				if (is.null(file)) { file <- "" }
				if(pretty == TRUE) {
					cat("____________________________\n", file=paste(file), append=TRUE)
					cat("Maximum likelihood estimation: rates:\n\n", file=paste(file), append=TRUE)
					cat("Lambda:          ", lambda, "\n", file=paste(file), append=TRUE)
					cat("Brownian variance (rate): ", mlparams$s2, "\n", file=paste(file), append=TRUE)
					cat("ML estimates of group relative rates :          ", max.lik.rate, "\n", file=paste(file), append=TRUE)
					cat("Lower confidence intervals for rates:", CIs.rate[,1], "\n", file=paste(file), append=TRUE)
					cat("Upper confidence intervals for rates:", CIs.rate[,2], "\n", file=paste(file), append=TRUE)
					cat("ML estimates of group means : ", mu, "\n\n", file=paste(file), append=TRUE)
					cat("Number of parameters: ", k, "\n", file=paste(file), append=TRUE)
					cat("Maximised log likelihood: ", max.lik, "\n", file=paste(file), append=TRUE)
					cat("  AIC = ", aic, "  \n", file=paste(file), append=TRUE)
					cat("  AICc = ", aicc, "  \n\n", file=paste(file), append=TRUE)
					cat("____________________________\n", file=paste(file), append=TRUE)
					cat("Comparison with single rate model\n", file=paste(file), append=TRUE)
					cat("Lambda (single rate model):          ", lambda_single, "\n", file=paste(file), append=TRUE)
					cat("Log likelihood (single rate): ", lik1, "\n", file=paste(file), append=TRUE)
					cat("LR statistic (ML rates vs single rate):", D, file=paste(file), append=TRUE)
					cat("  P = ", pval , file=paste(file), append=TRUE)
					cat(" df = ", k2, " \n", file=paste(file), append=TRUE)
					cat("  Single rate AIC = ", aic.rate1, "  \n", file=paste(file), append=TRUE)
					cat("  Single rate AICc = ", aicc.rate1, "  \n", file=paste(file), append=TRUE)
					cat("____________________________\n", file=paste(file), append=TRUE)
					}
				if(pretty == FALSE) {
					max.lik.rate.list <- list(BRvar=mlparams$s2, MLRate = max.lik.rate, LCI = CIs.rate[,1], UCI = CIs.rate[,2],  means=mu, nParam = k, lambda=lambda, Max.lik = max.lik, AIC = aic, AICc=aicc, LambdaSingle = lambda_single, Lik1 = lik1, LR = D, P = pval, df = k2, AIC.rate1=aic.rate1, AICc.rate1=aicc.rate1)
					return(max.lik.rate.list) }	
	}


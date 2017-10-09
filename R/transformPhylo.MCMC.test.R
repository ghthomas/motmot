#' @title Bayesian MCMC for models of trait evoluion (test version)
#' @description Fits Bayesian models for various models of continuous character evolution using a Markov Chain Monte Carlo (MCMC) approach
#' @param y A matrix of trait values.
#' @param phy An object of class "phylo" (see ape package).
#' @param model The model of trait evolution (see details).
#' @param mcmc.iteration Integer - the number of generations for which to run the MCMC chain
#' @param burn.in The proportion of the chain (as given by mcmc.iteration) which to discard as 'burn-in'
#' @param lowerBound Minimum value for parameter estimates
#' @param upperBound Maximum value for parameter estimates
#' @export


transformPhylo.MCMC.test <- function(y, phy, model, mcmc.iteration=1000, burn.in=0.1, lowerBound = NULL, upperBound = NULL, use.ML=FALSE, opt.accept.rate=TRUE, acceptance.sd=NULL) {
	
	bounds <- matrix(c(1e-08, 1, 1e-08, 1, 1e-08, 5, 1e-08, 20, 0, 1, 1e-08, 10, 1e-10, 100000), 7, 2, byrow = TRUE)
	rownames(bounds) <- c("kappa", "lambda", "delta", "alpha", "psi", "rate", "acdcRate")

	rtnorm <- function(n, mean, sd, a = -Inf, b = Inf){
		qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
	}

	if(model == "bm") {
			
			if(use.ML) {
				in.par <- transformPhylo.ML(y, phy, model="bm")
				input.value <- c(in.par[[1]])
			} else {
				input.value <- runif(1, 0, 1)
			}
				
			if (is.null(lowerBound)) {
				lowerBound <- bounds["rate", 1]
				}
			if (is.null(upperBound)) {
				upperBound <- bounds["rate", 2]
				}
						
			lik.model <- function(param) {
				return(likTraitPhylo.mcmc(y, phy, rate=param))
			}
			
			prior.uniform <- function(pram) {
				rate.prior <- dunif(pram[1], lowerBound[1], upperBound[1])
				return(rate.prior)
				}
				
			if(is.null(acceptance.sd)) {
				stn.dev <- 1
			} else {
				stn.dev <- acceptance.sd
			}
			
			name.param <- "brownianVariance"
		}
		
		
		if(model == "lambda") {
			
			if(use.ML) {
				in.par <- invisible(transformPhylo.ML(y, phy, model="lambda"))
				input.value <- c(in.par[[3]], in.par[[2]][1])
			} else {
				input.value <- runif(2, 0, 1)
			}
		
			if (is.null(lowerBound)) {
				lowerBound <- bounds["rate", 1]
				lowerBound[2] <- bounds["lambda", 1]
				}
			if (is.null(upperBound)) {
				upperBound <- bounds["rate", 2]
				upperBound[2] <-  bounds["lambda", 2]
				}
						
			lik.model <- function(pram) {
				lambda.phy <- transformPhylo(phy, model="lambda", y=y, lambda=pram[2])
				return(likTraitPhylo.mcmc(y, lambda.phy, rate=pram[1]))
			}
			
			prior.uniform <- function(pram) {
				rate.prior <- dunif(pram[1], lowerBound[1], upperBound[1])
				lambda.prior <- dunif(pram[2], lowerBound[2], upperBound[2])
				return(sum(rate.prior, lambda.prior))
				}
				
			if(is.null(acceptance.sd)) {
				stn.dev <- c(1, 0.001)
			} else {
				stn.dev <- acceptance.sd
			}	
				
			name.param <- c("brownianVariance", "Lambda")
		}
			
		
		if(model == "delta") {
			
			if(use.ML) {
				in.par <- invisible(transformPhylo.ML(y, phy, model="delta"))
				input.value <- c(in.par[[3]], in.par[[2]][1])
			} else {
				input.value <- runif(2, 0, 1)
			}
			
			if (is.null(lowerBound)) {
				lowerBound <- bounds["rate", 1]
				lowerBound[2] <- bounds["delta", 1]
				}
			if (is.null(upperBound)) {
				upperBound <- bounds["rate", 2]
				upperBound[2] <-  bounds["delta", 2]
				}
						
			lik.model <- function(pram) {
				delta.phy <- transformPhylo(phy, model="delta", y=y, delta=pram[2])
				return(likTraitPhylo.mcmc(y, delta.phy, rate=pram[1]))
			}
			
			prior.uniform <- function(pram) {
				rate.prior <- dunif(pram[1], lowerBound[1], upperBound[1])
				delta.prior <- dunif(pram[2], lowerBound[2], upperBound[2])
				return(sum(rate.prior, delta.prior))
				}
				
			if(is.null(acceptance.sd)) {
				stn.dev <- c(1, 0.001)
			} else {
				stn.dev <- acceptance.sd
			}	
			
			name.param <- c("brownianVariance", "Delta")
		}
		
		if(model == "kappa") {
			
			if(use.ML) {
				in.par <- invisible(transformPhylo.ML(y, phy, model="kappa"))
				input.value <- c(in.par[[3]], in.par[[2]][1])
			} else {
				input.value <- runif(2, 0, 1)
			}
			
		if (is.null(lowerBound)) {
				lowerBound <- bounds["rate", 1]
				lowerBound[2] <- bounds["kappa", 1]
				}
			if (is.null(upperBound)) {
				upperBound <- bounds["rate", 2]
				upperBound[2] <-  bounds["kappa", 2]
				}
						
			lik.model <- function(pram) {
				kappa.phy <- transformPhylo(phy, model="kappa", y=y, kappa=pram[2])
				return(likTraitPhylo.mcmc(y, kappa.phy, rate=pram[1]))
			}
			
			prior.uniform <- function(pram) {
				rate.prior <- dunif(pram[1], lowerBound[1], upperBound[1])
				kappa.prior <- dunif(pram[2], lowerBound[2], upperBound[2])
				return(sum(rate.prior, kappa.prior))
				}
				
			if(is.null(acceptance.sd)) {
				stn.dev <- c(1, 0.001)
			} else {
				stn.dev <- acceptance.sd
			}	

			name.param <- c("brownianVariance", "Kappa")
		}
		
		if(model == "OU") {
			
			if(use.ML) {
				in.par <- invisible(transformPhylo.ML(y, phy, model="OU", modelCIs=F))
				input.value <- c(in.par[[3]], in.par[[2]][1])
			} else {
				input.value <- runif(2, 0, 1)
			}
			
			if (is.null(lowerBound)) {
				lowerBound <- bounds["rate", 1]
				lowerBound[2] <- bounds["alpha", 1]
				}
			if (is.null(upperBound)) {
				upperBound <- bounds["rate", 2]
				upperBound[2] <-  bounds["alpha", 2]
				}
						
			lik.model <- function(pram) {
				ou.phy <- transformPhylo(phy, model="OU", y=y, alpha=pram[2])
				return(likTraitPhylo.mcmc(y, ou.phy, rate=pram[1]))
			}
			
			prior.uniform <- function(pram) {
				rate.prior <- dunif(pram[1], lowerBound[1], upperBound[1])
				ou.prior <- dunif(pram[2], lowerBound[2], upperBound[2])
				return(sum(rate.prior, ou.prior))
				}
				
			if(is.null(acceptance.sd)) {
				stn.dev <- c(1, 0.001)
			} else {
				stn.dev <- acceptance.sd
			}	
				
			name.param <- c("brownianVariance", "OU")
		}
		
		if(model == "ACDC") {
			
			if(use.ML) {
				in.par <- invisible(transformPhylo.ML(y, phy, model="ACDC"))
				input.value <- c(in.par[[3]], in.par[[2]][1])
			} else {
				input.value <- runif(2, 0, 1)
			}

			rootBranchingTime <- nodeTimes(phy)[1,1]
			if (is.null(lowerBound)) {
				lowerBound <- bounds["rate", 1]
				lowerBound[2] <- log(bounds["acdcRate", 1]) / rootBranchingTime
				}
			if (is.null(upperBound)) {
				upperBound <- bounds["rate", 2]
				upperBound[2] <- log(bounds["acdcRate", 2]) / rootBranchingTime
				}
						
			lik.model <- function(pram) {
				acdc.phy <- transformPhylo(phy, model="ACDC", y=y, acdcRate=pram[2])
				return(likTraitPhylo.mcmc(y, acdc.phy, rate=pram[1]))
			}
			
			prior.uniform <- function(pram) {
				rate.prior <- dunif(pram[1], lowerBound[1], upperBound[1])
				acdc.prior <- dunif(pram[2], lowerBound[2], upperBound[2])
				return(sum(rate.prior, acdc.prior))
				}
				
			if(is.null(acceptance.sd)) {
				stn.dev <- c(1, 1)
			} else {
				stn.dev <- acceptance.sd
			}		

			name.param <- c("brownianVariance", "ACDC.rate")
		}
		
		if(model == "psi") {
			
			if(use.ML) {
				in.par <- invisible(transformPhylo.ML(y, phy, model="psi"))
				input.value <- c(in.par[[3]], in.par[[2]][1])
			} else {
				input.value <- runif(2, 0, 1)
			}

			if (is.null(lowerBound)) {
				lowerBound <- bounds["rate", 1]
				lowerBound[2] <- bounds["psi", 1]
				}
			if (is.null(upperBound)) {
				upperBound <- bounds["rate", 2]
				upperBound[2] <- bounds["psi", 2]
				}
						
			lik.model <- function(pram) {
				psi.phy <- transformPhylo(phy, model="psi", y=y, psi=pram[2])
				return(likTraitPhylo.mcmc(y, psi.phy, rate=pram[1]))
			}
			
			prior.uniform <- function(pram) {
				rate.prior <- dunif(pram[1], lowerBound[1], upperBound[1])
				psi.prior <- dunif(pram[2], lowerBound[2], upperBound[2])
				return(sum(rate.prior, psi.prior))
				}
				
			if(is.null(acceptance.sd)) {
				stn.dev <- c(1, 0.1)
			} else {
				stn.dev <- acceptance.sd
			}	
			name.param <- c("brownianVariance", "psi")
		}
		
		model.posterior <- function(pram) return (lik.model(pram) + prior.uniform(pram))		
		
		motmot.mcmc <- function(input.value, iterations, stn.dev, silent=FALSE)  {
		
				propose.mcmc <- function(pram) {
					return(sapply(1:length(stn.dev) , function(x) rtnorm(1, pram[x], sd=stn.dev[x], lowerBound[x], upperBound[x])))
				}
					
			mcmc.chain <- matrix(input.value, nrow=1)
    				for (i in 1:iterations) {
        				proposed.move <- propose.mcmc(mcmc.chain[i,])
        				chain.prob <- exp(model.posterior(proposed.move) - model.posterior(mcmc.chain[i,]))
        				if (runif(1) < chain.prob) {
        					mcmc.chain <- rbind(mcmc.chain, proposed.move)
        				} else {
        					mcmc.chain <- rbind(mcmc.chain, mcmc.chain[i,])
        					}
					if(!silent) {
						cat("\r", "MCMC progress:", sprintf("%.4f", i/iterations * 100), "%")
						}
    					}
    				return(mcmc.chain)
				}
				
		if(opt.accept.rate) {
			cat("optimising acceptance proposal")
			cat("\n", " ")	
			cat("\n")	
			count <- 1
			stn.dev.2 <- stn.dev
			old.ratio <- 10		
		
			opt.done <- FALSE
			while(opt.done == FALSE) {
				chain <- motmot.mcmc(input.value, mcmc.iteration / 10, stn.dev.2, silent=TRUE)
				burnIn <- ceiling((mcmc.iteration / 10) * burn.in)
				acceptance <- 1 - mean(duplicated(chain[-(1:burnIn),]))		
				diff.to.accept <- acceptance - 0.3
				new.ratio <- abs(acceptance - 0.3)
				if(new.ratio < old.ratio) {
					stn.dev <- stn.dev.2
					old.ratio <- new.ratio
					cat("\r", count, " better acceptance found")
					count <- count + 1
					}					
				if(abs(acceptance - 0.3) < 0.02) {
					opt.done <- TRUE
					cat(" finished.")
					}
				stn.dev.2 <- sapply(1:length(stn.dev), function(x) rtnorm(1, stn.dev[x], 0.1, 1e-8, Inf))
				}		
			}			
		cat("\n", " ")		
		chain <- motmot.mcmc(input.value, mcmc.iteration, stn.dev=stn.dev)
		burnIn <- ceiling(mcmc.iteration * burn.in)
		acceptance.1 <- 1 - mean(duplicated(chain[-(1:burnIn), 1]))
		post.burn.in <- chain[-c(1:burnIn), ]
		
		if(dim(chain)[2] == 1) {
			ess.mcmc <- effectiveSize(post.burn.in)
			median.mcmc <- median(post.burn.in)
			hpd.mcmc <- quantile(post.burn.in, c(0.025, 0.975))
			names(ess.mcmc) <- names(median.mcmc) <- name.param
			names(hpd.mcmc) <- c("lower 95% HPD", "upper 95% HPD")
		} else {
			ess.mcmc <- apply(post.burn.in, 2, effectiveSize)
			median.mcmc <- apply(chain[-c(1:burnIn), ], 2, median)
			hpd.mcmc <- apply(chain[-c(1:burnIn), ], 2, function(x) quantile(x, c(0.025, 0.975)))
			names(ess.mcmc) <- names(median.mcmc) <- name.param
			colnames(hpd.mcmc) <- name.param
		}
		
		cat("\n")
		output.mcmc <- list(median.mcmc, hpd.mcmc, ess.mcmc, acceptance.1, chain)
		names(output.mcmc) <- c("median", "95.HPD", "ESS", "acceptance.ratio", "mcmc.chain")
		print(output.mcmc[1:4])
		invisible(return(output.mcmc))
	}
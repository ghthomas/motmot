likRatePhylo <-
function(rateData, rate=NULL, common.mean=FALSE, lambda.est=TRUE, lambda=1, meserr=FALSE, sigmaScale=NULL) {
	
		if(is.null(rate))  { rate <- c(rep(1,length(rateData$Vmat))) } else { rate <- rate }	

			y <- rateData$y
			x <- as.factor(rateData$x)
			if (common.mean==FALSE) {k <- nlevels(x)} else {k <- 1}
	
			V <- transformRateMatrix(rateData, rate)

			x <- make.anc(y, x)
						
			if (!lambda.est & !meserr) {
				logDetV <- determinant(V)$modulus
				mu <- phyloMean(rateData, rate, common.mean, lambda.est, lambda)
				s2 <- phyloVar(rateData, rate, common.mean, lambda.est, lambda)
				n <- length(x[,1])
				ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(s2) - logDetV / 2.0 - (n - k)/2.0
				lik.RatePhylo <- ( list(ll = ll, mu = mu, s2 = s2) )
			}

			if (lambda.est & !meserr) {
				v.temp <- V
				diag(v.temp) <- rep(0, dim(V)[1])
				V.lam <- lambda*v.temp
				diag(V.lam) <- diag(V)
				V <- V.lam
				logDetV <- determinant(V)$modulus
				mu <- phyloMean(rateData, rate, common.mean, lambda.est, lambda)
				s2 <- phyloVar(rateData, rate, common.mean, lambda.est, lambda)
				n <- length(x[,1])
				ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(s2) - logDetV / 2.0 - (n - k)/2.0
				lik.RatePhylo <- ( list(ll = ll, mu = mu, s2 = s2) )
			}

	
			if (!lambda.est & meserr) {
				
				if (is.null(sigmaScale)) { stop("Estimate of sigmaScale required") }
				
				if(common.mean==FALSE) {x <- x} else { x <- rep(1, length(x[,1]))}

				diag(V) <- diag(V) + rateData$meserr/sigmaScale
				logDetV <- determinant(V)$modulus
				
				iV <- solve(V)
				xVix <- crossprod(x, iV %*% x)
				xViy <- crossprod(x, iV %*% y)
				mu <- solve(xVix) %*% xViy 

				e <- y - x %*% mu
				s2 <- crossprod(e, iV %*% e)
				n <- length(y) 
				phylo.var <- ( s2 / (n - k) )
				
				
				ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(phylo.var) - logDetV / 2.0 - (n - k)/2.0
				lik.RatePhylo <- ( list(ll = ll, mu = mu, s2 = phylo.var) )
			}
	
	
			if (lambda.est & meserr) {
		
				if (is.null(sigmaScale)) { stop("Estimate of sigmaScale required") }
		
				if(common.mean==FALSE) {x <- x} else { x <- rep(1, length(x[,1]))}
		
				
				v.temp <- V
				diag(v.temp) <- rep(0, dim(V)[1])
				V.lam <- lambda*v.temp
				diag(V.lam) <- diag(V)
				V <- V.lam
				diag(V) <- diag(V) + rateData$meserr/sigmaScale
				logDetV <- determinant(V)$modulus
		
				iV <- solve(V)
				xVix <- crossprod(x, iV %*% x)
				xViy <- crossprod(x, iV %*% y)
				mu <- solve(xVix) %*% xViy 
		
				e <- y - x %*% mu
				s2 <- crossprod(e, iV %*% e)
				n <- length(y) 
				phylo.var <- ( s2 / (n - k) )
		
		
				ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(phylo.var) - logDetV / 2.0 - (n - k)/2.0
				lik.RatePhylo <- ( list(ll = ll, mu = mu, s2 = phylo.var, lambda=lambda) )
			}
	
	
				return(lik.RatePhylo)
	}







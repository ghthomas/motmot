make.likRatePhylo <-
function(rateData, fixed, common.mean=FALSE, lambda.est, meserr) {
	
	op <- c(fixed, !lambda.est, !meserr)
	
	function(params){
		
	
		op[c(!fixed, lambda.est, meserr)] <- params[c(!fixed, lambda.est, meserr)]
		

		y <- rateData$y
		x <- as.factor(rateData$x)
		if (common.mean==FALSE) {k <- nlevels(x)} else {k <- 1}
		
		V <- rateData$Vmat
		v1 <- V[[1]]
		nV <- length(rateData$Vmat)
		
		rateMats <- vector(mode="list", length = nV)
		retMat <- matrix(0, nrow = dim(v1)[1], ncol = dim(v1)[2])
		
		for(i in 1:nV) {
			rateMats[[i]] <- op[i] * V[[i]]  
			retMat <- retMat + rateMats[[i]]
		}
										
			x <- make.anc(y, x)
			
			if(common.mean==FALSE) {
				x <- x} else { x <- rep(1, length(x[,1]))}
									
				V <- retMat
				v.temp <- V
				diag(v.temp) <- rep(0, dim(V)[1])
			
				V.lam <- op[nV+1]*v.temp
				diag(V.lam) <- diag(V)
				V <- V.lam
				if (meserr==TRUE) {	diag(V) <- diag(V) + rateData$meserr/op[nV+2] } 

				logDetV <- determinant(V)$modulus
				
				iV <- solve(V)
				xVix <- crossprod(x, iV %*% x)
				xViy <- crossprod(x, iV %*% y)
				mu <- solve(xVix) %*% xViy 
				
				e <- y - x %*% mu
				s2 <- crossprod(e, iV %*% e)
				n <- length(y) 
				phylo.var <- ( s2 / (n - k) )
				
				n <- length(y)
				ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(phylo.var) - logDetV / 2.0 - (n - k)/2.0
				
				lik.RatePhylo <- ( list(ll = ll, mu = mu, phylo.var = phylo.var, lambda=op[nV+2]) )
				return(-1 * lik.RatePhylo$ll)	
		}
	
	}


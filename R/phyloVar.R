phyloVar <-
function(rateData, rate=NULL, common.mean=FALSE, lambda.est=TRUE, lambda=1, meserr=FALSE) {
	
		if(is.null(rate))  { rate <- c(rep(1,length(rateData$Vmat))) } else { rate <- rate }	
	
		if(length(rate) != length(rateData$Vmat)){stop("The number of rates defined differs from the number of rate matrices")}
			
		y <- rateData$y
		x <- as.factor(rateData$x)
						
		if (common.mean==FALSE) {k <- nlevels(x)} else {k <- 1}
			
		V <- transformRateMatrix(rateData, rate)
		if (lambda.est & !meserr) {
			v.temp <- V
			diag(v.temp) <- rep(0, dim(V)[1])
			V.lam <- lambda*v.temp
			diag(V.lam) <- diag(V)
			V <- V.lam
		}
		
		x <- make.anc(y, x)
			
		if(common.mean==FALSE) {x <- x} else { x <- rep(1, length(x[,1]))}
		
		mu <- phyloMean(rateData, rate, common.mean=common.mean, lambda.est, lambda)
			
			iV <- solve(V)
			e <- y - x %*% mu
			s2 <- crossprod(e, iV %*% e)
			n <- length(y) 
			phylo.var <- ( s2 / (n - k) )
			return(phylo.var)
			}


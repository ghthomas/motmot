#' @title Calculation of phylogenetically corrected mean.
#' @description This function calculates the phylogenetic mean of the data given the tree and model of evolution
#' @param rateData an object of class "rateData"
#' @param rate a vector of relative rate parameters. The length of the vector is equal to the number of rates being estimated. If \code{rate=NULL} then rates are equal.
#' @param common.mean a logical specififying whether each rate category should have its own mean (\code{common.mean=FALSE}) or all categories should have the same mean (\code{common.mean=FALSE}). See Thomas et al. (2009) for a discussion on the impact of assumptions about mean on rate estimates.
#' @param lambda.est Logical. Fit Pagel's lambda.
#' @param lambda Numeric value for lambda from 0-1.
#' @param meserr Logical. Include measurement error.
#' @return mu phylogenetically corrected mean
#' @references Thomas GH, Freckleton RP, & Szekely T. 2006. Comparative analyses of the influence of developmental mode on phenotypic diversification rates in shorebirds. Proceedings of the Royal Society B 273, 1619-1624.
#' Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments and island effects. Evolution 63, 2017-2030.
#' @author Gavin Thomas, Rob Freckleton
#' @note The means are output as treatment contrasts.
#' @export

phyloMean <-
function(rateData, rate=NULL, common.mean=FALSE, lambda.est=TRUE, lambda=1, meserr=FALSE) {
	
		if(is.null(rate))  { rate <- c(rep(1,length(rateData$Vmat))) } else { rate <- rate }	
				
		if(length(rate) != length(rateData$Vmat)){stop("The number of rates defined differs from the number of rate matrices")}
	
		y <- rateData$y
		x <- as.factor(rateData$x)

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

			iV <- solve(V)
			xVix <- crossprod(x, iV %*% x)
			xViy <- crossprod(x, iV %*% y)
			mu <- solve(xVix) %*% xViy 
			return(mu)
		}


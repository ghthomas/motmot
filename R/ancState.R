#' @export

ancState <- function (phy, y) {

			V <- vcv(phy)
			x <- rep(1, Ntip(phy))
			iV <- solve(V)
			xVix <- crossprod(x, iV %*% x)
			xViy <- crossprod(x, iV %*% y)
			mu <- solve(xVix) %*% xViy
			return(as.numeric(mu))
  }  

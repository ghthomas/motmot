#' @export

likTraitPhylo.mcmc <- function (y, phy, rate, meserr = NULL) {

    n <- length(phy$tip.label)
    k <- ncol(y)
    y <- as.matrix(y[phy$tip.label, ])
    contrasts <- pic.motmot(y, phy)
    rawVariances <- c(contrasts$contr[, 2], contrasts$V)
    rawContrasts <- matrix(c(contrasts$contr[, 1], 0), ncol=1)
    brCov <- as.matrix(rate)
    iW <- solve(brCov)
	addCons <- 0
    for (i in 1:n) {
        ui <-rawContrasts[i, ]
        cpUI <- iW %*% ui
        addCons <- addCons  + ( ui %*% cpUI ) / rawVariances[i]
	}	
	
    return(-0.5 * (n * k * log(2 * pi) + n * log(det(brCov)) + k * sum(log(rawVariances)) + addCons))
}

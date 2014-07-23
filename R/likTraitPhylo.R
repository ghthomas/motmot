likTraitPhylo<-function (y, phy, meserr = NULL, covPIC = TRUE)
{
    if (is.matrix(y) == FALSE) {
        stop("Trait data must be a matrix with taxon names as row names")
    }
    n <- length(phy$tip.label)
    k <- ncol(y)
    phy <- reorder(phy, order = "pruningwise")
    y <- as.matrix(y[phy$tip.label, ])
    contrasts <- apply(y, 2, pic.motmot, phy = phy)
    rawVariances <- c(contrasts[[1]]$contr[, 2], contrasts[[1]]$V)
    rawContrasts <- matrix(NA, nrow = n, ncol = ncol(y))
    for (i in 1:k) {
        rawContrasts[, i] <- c(contrasts[[i]]$contr[, 1], 0)
    }
    brCov <- matrix(NA, nrow = ncol(y), ncol = ncol(y))
    for (i in 1:k) {
        for (j in 1:k) {
			if(i==j|covPIC==TRUE){
                brCov[j, i] <- brCov[i, j] <- crossprod(rawContrasts[,
                j]/sqrt(rawVariances), rawContrasts[, i]/sqrt(rawVariances))/(n - 1)
			}else{
                brCov[j,i]<-0
			}
        }
    }
    iW <- solve(brCov)
    addCons <- 0
    for (i in 1:n) {
        ui <- matrix(rawContrasts[i, ])
        addCons <- addCons + crossprod(ui, iW %*% ui)/rawVariances[i]
    }
    logLikelihood <- -0.5 * (n * k * log(2 * pi) + n * log(det(brCov)) +
    k * sum(log(rawVariances)) + addCons)
    return(list(brownianVariance = brCov, logLikelihood = logLikelihood))
}

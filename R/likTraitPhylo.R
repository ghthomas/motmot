likTraitPhylo <- function (y, phy, meserr = NULL) 
{
    if (is.matrix(y) == FALSE) {
        stop("Trait data must be a matrix with taxon names as row names")
    }
    n <- length(phy$tip.label)
    k <- ncol(y)
#phy <- reorder(phy, order = "pruningwise")
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
            brCov[j, i] <- brCov[i, j] <- crossprod(rawContrasts[, 
													j]/sqrt(rawVariances), rawContrasts[, i]/sqrt(rawVariances))/(n - 
																												  1)
        }
    }
    iW <- solve(brCov)
	addCons <- 0
    
    for (i in 1:n) {
        ui <-rawContrasts[i, ]
        cpUI <- iW %*% ui
        addCons <- addCons  + ( ui %*% cpUI ) / rawVariances[i]
	}	
	
#	allConstrasts <- cbind(rawContrasts, rawVariances)
#   allConstrasts <- matrix(allConstrasts, nrow=length(rawVariances))
#	allConstrastsList <- list()
#	for (i in 1:length(allConstrasts[,1])) {allConstrastsList[[i]] <- allConstrasts[i,]}
	
#	foo <- function(x, iW, k) {
#		ui <- x[1:length(x)-1]
#		cpUI <- iW %*% ui
#		return(( ui %*% cpUI ) / x[length(x)])
#		}
	
#	addCons <- sum(unlist(mclapply(allConstrastsList, FUN=foo, iW=iW, k=k, mc.cores=n.cores)))
	
	
    logLikelihood <- -0.5 * (n * k * log(2 * pi) + n * log(det(brCov)) + 
							 k * sum(log(rawVariances)) + addCons)
    return(list(brownianVariance = brCov, logLikelihood = logLikelihood))
}

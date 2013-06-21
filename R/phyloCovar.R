phyloCovar <- function(x, phy, estimator="unbiased") {
		
		if (is.matrix(x)==FALSE) { stop("Trait data must be a matrix with taxon names as row names")}
		
		n <- length(phy$tip.label)
		phy <- reorder(phy, order = "pruningwise")
	
		x <- as.matrix(x[phy$tip.label,])
		
		contrasts <- apply(x, MARGIN=2, FUN=pic.motmot, phy=phy)
		rawVariances <- c(contrasts[[1]]$contr[,2], contrasts[[1]]$V)
		rawContrasts <- matrix(NA, nrow=n, ncol=ncol(x))
		
		for (i in 1:ncol(x)) {
			rawContrasts[,i] <- c(contrasts[[i]]$contr[,1],0)
			}
			
		brCov <- matrix(NA, nrow=ncol(x), ncol=ncol(x))
		
	if (estimator=="unbiased") {
		for (i in 1:ncol(x)) {
			for (k in 1:ncol(x)) {	
				brCov[k,i] <- brCov[i,k] <- crossprod(rawContrasts[,k]/sqrt(rawVariances), rawContrasts[,i]/sqrt(rawVariances)) / (n-1)
								}}
	}
	
	if (estimator=="ML") {
		for (i in 1:ncol(x)) {
			for (k in 1:ncol(x)) {	
				brCov[k,i] <- brCov[i,k] <- crossprod(rawContrasts[,k]/sqrt(rawVariances), rawContrasts[,i]/sqrt(rawVariances)) / (n)
			}}
	}
		
		return(brCov)
		}
		


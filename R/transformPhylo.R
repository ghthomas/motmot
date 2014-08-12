transformPhylo<-function (phy, model = NULL, meserr = NULL, y = NULL, kappa = NULL,  lambda = NULL, delta = NULL, alpha = NULL, psi = NULL, la = NULL, nodeIDs = NULL, rateType = NULL, branchRates = NULL, cladeRates = NULL, branchLabels = NULL)
{
    if (is.null(meserr) == FALSE) {
        if (dim(y)[2] > 1) {
            meserr <- NULL
            (stop("Measurement error can only be included for univariate models. Set meserr to NULL."))
        }
    }
    n <- length(phy$tip.label)
    switch(model, bm = {
        if (is.null(meserr) == FALSE) {
            height <- max(branching.times(phy))
            interns <- which(phy$edge[, 2] > n)
            externs <- which(phy$edge[, 2] <= n)
            phy$edge.length[externs] <- phy$edge.length[externs] +
            (meserr^2)/(var(y)/height)[1]
        } else {
            phy <- phy
        }
    }, kappa = {
        if (is.null(meserr) == FALSE) {
            height <- max(branching.times(phy))
            interns <- which(phy$edge[, 2] > n)
            externs <- which(phy$edge[, 2] <= n)
        }
        phy$edge.length <- phy$edge.length^kappa
        if (is.null(meserr) == FALSE) {
            phy$edge.length[externs] <- phy$edge.length[externs] +
            (meserr^2)/(var(y)/height)[1]
        }
    }, lambda = {
        if (is.null(meserr) == FALSE) {
            height <- max(branching.times(phy))
            interns <- which(phy$edge[, 2] > n)
            externs <- which(phy$edge[, 2] <= n)
        }
        if (is.ultrametric(phy)) {
            rootOrig <- max(branching.times(phy))
            tips <- match(c(1:Ntip(phy)), phy$edge[, 2])
            phy$edge.length <- phy$edge.length * lambda
            phy$edge.length[tips] <- phy$edge.length[tips] +
            (rootOrig * (1 - lambda))
        }
        if (is.ultrametric(phy) == FALSE) {
            tips <- match(c(1:Ntip(phy)), phy$edge[, 2])
            cladeMat <- clade.matrix(phy)
            branchHeights <- rep(NA, Ntip(phy))
            for (i in 1:Ntip(phy)) {
                branchHeights[i] <- sum(cladeMat$edge.length[cladeMat$clade.matrix[,
                i] == 1])
            }
            phy$edge.length <- phy$edge.length * lambda
            phy$edge.length[tips] <- phy$edge.length[tips] +
            (branchHeights * (1 - lambda))
        }
        if (is.null(meserr) == FALSE) {
            phy$edge.length[externs] <- phy$edge.length[externs] +
            (meserr^2)/(var(y)/height)[1]
        }
    }, delta = {
        if (is.null(meserr) == FALSE) {
            height <- max(branching.times(phy))
            interns <- which(phy$edge[, 2] > n)
            externs <- which(phy$edge[, 2] <= n)
        }
        times <- branching.times(phy)
        times <- max(times) - times
        tips <- length(phy$tip.label)
        res <- phy
        for (i in 1:length(phy$edge.length)) {
            bl <- phy$edge.length[i]
            age <- times[phy$edge[i, 1] - tips]
            res$edge.length[i] <- (age + bl)^delta - age^delta
        }
        phy <- res
        if (is.null(meserr) == FALSE) {
            phy$edge.length[externs] <- phy$edge.length[externs] +
            (meserr^2)/(var(y)/height)[1]
        }
    }, free = {
        if (is.null(meserr) == FALSE) {
            height <- max(branching.times(phy))
            interns <- which(phy$edge[, 2] > n)
            externs <- which(phy$edge[, 2] <= n)
        }
        branchRates <- branchRates + (1 - min(branchRates))
        phy$edge.length <- phy$edge.length * branchRates
        if (is.null(meserr) == FALSE) {
            phy$edge.length[externs] <- phy$edge.length[externs] +
            (meserr^2)/(var(y)/height)[1]
        }
    }, clade = {
        if (is.null(meserr) == FALSE) {
            height <- max(branching.times(phy))
            interns <- which(phy$edge[, 2] > n)
            externs <- which(phy$edge[, 2] <= n)
        }
        if (is.null(rateType)) {
            rateType <- rep("clade", length(nodeIDs))
        } else {
            rateType <- rateType
        }
        branchShiftNms <- branchShiftID <- cladeShiftNms <- cladeShiftID <- NULL
        cladeMembers <- matrix(0, ncol = length(nodeIDs), nrow = length(phy$edge[,
        1]))
        shiftType <- data.frame(rateType, nodeIDs, cladeRates)
        colnms <- paste(shiftType[, 1], shiftType[, 2], sep = "")
        shiftType <- data.frame(shiftType, colnms)
        if (sum(shiftType[, 1] == "clade") > 0) {
            cladeShiftID <- shiftType[shiftType[, 1] == "clade",
            2]
            cladeShiftNms <- as.character(shiftType[shiftType[,
            1] == "clade", 4])
            cladeMembers[, 1:length(cladeShiftID)] <- cladeIdentity(phy = phy,
            nodeIDs = cladeShiftID)
        }
        if (sum(shiftType[, 1] == "branch") > 0) {
            branchShiftID <- shiftType[shiftType[, 1] == "branch",
            2]
            branchShiftNms <- as.character(shiftType[shiftType[,
            1] == "branch", 4])
        }
        if (is.null(branchShiftNms) == FALSE) {
            colnames(cladeMembers) <- c(cladeShiftNms, rep(NA,
            length(branchShiftNms)))
        }
        if (is.null(branchShiftNms) == TRUE) {
            colnames(cladeMembers) <- cladeShiftNms
        }
        if (sum(shiftType[, 1] == "branch") > 0) {
            for (i in 1:length(branchShiftID)) {
                branchID <- which(phy$edge[, 2] == branchShiftID[i])
                cladeMembers[branchID, ] <- 0
                cladeMembers[branchID, length(cladeShiftID) +
                i] <- 1
                colnames(cladeMembers)[length(cladeShiftID) +
                i] <- branchShiftNms[i]
            }
        }
        cladeMembers <- as.matrix(cladeMembers[, match(shiftType[,
        4], colnames(cladeMembers))])
        for (i in 1:length(cladeRates)) {
            phy$edge.length[cladeMembers[, i] == 1] <- phy$edge.length[cladeMembers[,
            i] == 1] * cladeRates[i]
        }
        if (is.null(meserr) == FALSE) {
            phy$edge.length[externs] <- phy$edge.length[externs] +
            (meserr^2)/(var(y)/height)[1]
        }
    }, OU = {
        if (is.ultrametric(phy)==FALSE) {stop("This implementation of the OU model can only be applied to ultrametric trees. See Slater, 2014 DOI: 10.1111/2041-210X.12201")
        if (is.null(meserr) == FALSE) {
            height <- max(branching.times(phy))
            interns <- which(phy$edge[, 2] > n)
            externs <- which(phy$edge[, 2] <= n)
        }
        times <- branching.times(phy)
        names(times) <- (Ntip(phy) + 1):(Ntip(phy) + Nnode(phy))
        Tmax <- times[1]
        age <- times[match(phy$edge[,1], names(times))]
        bl <- phy$edge.length
        t1 <- max(times) - age
        t2 <- t1 + bl
        
        phy$edge.length <- (1/(2 * alpha)) * exp(-2 * alpha * (Tmax - t2)) * (1 - exp(-2 * alpha * t2)) - (1/(2 * alpha)) * exp(-2 * alpha * (Tmax - t1)) * (1 - exp(-2 * alpha * t1))
        
        if (is.null(meserr) == FALSE) {
            phy$edge.length[externs] <- phy$edge.length[externs] +
            (meserr^2)/(var(y)/height)[1]
        }
    }, psi = {
        if (is.null(meserr) == FALSE) {
            height <- max(branching.times(phy))
            interns <- which(phy$edge[, 2] > n)
            externs <- which(phy$edge[, 2] <= n)
        }
		if(is.null(phy$Sobs)){
			Stot<-rep(1, length(phy$edge.length))
		}else{
			Stot<-phy$Sobs
		}
		if(!is.null(phy$Shid)) Stot<-Stot+phy$Shid
        phy2 <- phy
        phy2$edge.length <- (psi/(2*la)) * Stot + (1 - psi) * phy$edge.length
        phy <- phy2
        if (is.null(meserr) == FALSE) {
            phy$edge.length[externs] <- phy$edge.length[externs] +
            (meserr^2)/(var(y)/height)[1]
        }
    }, multipsi = {
        if (is.null(meserr) == FALSE) {
            height <- max(branching.times(phy))
            interns <- which(phy$edge[, 2] > n)
            externs <- which(phy$edge[, 2] <= n)
        }
		if(is.null(phy$Sobs)){
			Stot<-rep(1, length(phy$edge.length))
		}else{
			Stot<-phy$Sobs
		}
		if(!is.null(phy$Shid)) Stot<-Stot+phy$Shid
        phy2 <- phy
		states<-levels(factor(branchLabels))
		for(i in 1:length(states)){
            phy2$edge.length[branchLabels==states[i]] <- (psi[i]/(2*la)) * Stot[branchLabels==states[i]] + (1 - psi[i]) * phy$edge.length[branchLabels==states[i]]
		}
        phy <- phy2
        if (is.null(meserr) == FALSE) {
            phy$edge.length[externs] <- phy$edge.length[externs] + (meserr^2)/(var(y)/height)[1]
        }
    })
    return(phy)
}
transformPhylo <- function (phy, model = NULL, y = NULL, kappa = NULL, 
lambda = NULL, delta = NULL, alpha = NULL, psi = NULL, nodeIDs = NULL, 
rateType = NULL, branchRates = NULL, cladeRates = NULL, cladeMembersObj = NULL, meserr=NULL, sigmaScale=1) 
{

    n <- length(phy$tip.label)
	
	if (is.null(meserr)) { meserr <- rep(0, Ntip(phy)) }
	if (length(meserr)==1) { meserr <- (meserr * y)^2 }
	if (length(meserr)>1) { meserr <- meserr^2 }

		
		
    switch(model, bm = {
		   if (is.null(meserr) == FALSE) {
		   height <- max(branching.times(phy))
		   interns <- which(phy$edge[, 2] > n)
		   externs <- which(phy$edge[, 2] <= n)
		   phy$edge.length[externs] <- phy$edge.length[externs] + (meserr^2)*sigmaScale
		   } else {
		   phy <- phy
		   }
		   
		   }, kappa = {
		   height <- max(branching.times(phy))
		   if (is.null(meserr) == FALSE) {
		   height <- max(branching.times(phy))
		   interns <- which(phy$edge[, 2] > n)
		   externs <- which(phy$edge[, 2] <= n)
		   }
		   phy$edge.length <- phy$edge.length^kappa
		   if (is.null(meserr) == FALSE) {
		   phy$edge.length[externs] <- phy$edge.length[externs] + (meserr^2)*sigmaScale
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
		   height <- max(branching.times(phy))
		   if (is.null(meserr) == FALSE) {
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
		   phy$edge.length[externs] <- phy$edge.length[externs] + (meserr^2)*sigmaScale
		   }
		   
		   phy$edge.length <- phy$edge.length * (height/ max(branching.times(phy)))
		   
		   
		   }, free = {
		   if (is.null(meserr) == FALSE) {
		   height <- max(branching.times(phy))
		   interns <- which(phy$edge[, 2] > n)
		   externs <- which(phy$edge[, 2] <= n)
		   }
		   branchRates <- branchRates + (1 - min(branchRates))
		   phy$edge.length <- phy$edge.length * branchRates
		   if (is.null(meserr) == FALSE) {
		   phy$edge.length[externs] <- phy$edge.length[externs] + (meserr^2)*sigmaScale
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
																   nodeIDs = cladeShiftID, cladeMembersObj = cladeMembersObj)
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
		   if (sum(meserr) != 0) {
		   phy$edge.length[externs] <- phy$edge.length[externs] + 
		   (meserr^2)/(var(y)/height)[1]
		   }
		   }, OU = {
		   
		   times <- branching.times(phy)
		   height <- max(times)
		   
		   if (is.null(meserr) == FALSE) {
		   interns <- which(phy$edge[, 2] > n)
		   externs <- which(phy$edge[, 2] <= n)
		   }
		   names(times) <- (Ntip(phy) + 1):(Ntip(phy) + Nnode(phy))
		   Tmax <- times[1]
		   phy2 <- phy
		   for (i in 1:length(phy$edge.length)) {
		   bl <- phy$edge.length[i]
		   age <- times[which(names(times) == phy$edge[i, 1])]
		   t1 <- max(times) - age
		   t2 <- t1 + bl
		   phy2$edge.length[i] <- (1/(2 * alpha)) * exp(-2 * 
														alpha * (Tmax - t2)) * (1 - exp(-2 * alpha * 
																						t2)) - (1/(2 * alpha)) * exp(-2 * alpha * (Tmax - 
																																   t1)) * (1 - exp(-2 * alpha * t1))
		   }
		   phy <- phy2
		   if (is.null(meserr) == FALSE) {
		   phy$edge.length[externs] <- phy$edge.length[externs] + 
		   (meserr^2)/(var(y)/height)[1]
		   }
		   
		   phy$edge.length <- phy$edge.length * (height/ max(branching.times(phy)))
		   
		   
		   }, psi = {
		   if (is.null(meserr) == FALSE) {
		   height <- max(branching.times(phy))
		   interns <- which(phy$edge[, 2] > n)
		   externs <- which(phy$edge[, 2] <= n)
		   }
		   bdrates <- birthdeath(phy)
		   mu <- bdrates$para["d/b"] * bdrates$para["b-d"]
		   lambda <- bdrates$para["b-d"] + mu
		   probNoDescendents <- mu * (exp((lambda - mu) * phy$edge.length) - 
									  1)/((lambda * exp((lambda - mu) * phy$edge.length)) - 
										  mu)
		   sh <- (lambda * probNoDescendents) * phy$edge.length
		   phy2 <- phy
		   phy2$edge.length <- (psi/lambda) * (phy$edge.length^0 + 
											   sh) + (1 - psi) * phy$edge.length
		   phy <- phy2
		   if (is.null(meserr) == FALSE) {
		   phy$edge.length[externs] <- phy$edge.length[externs] + 
		   (meserr^2)/(var(y)/height)[1]
		   }
		   })
    return(phy)
}

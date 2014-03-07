transformPhylo.ML <- function (y, phy, model = NULL, modelCIs = TRUE, 
nodeIDs = NULL, rateType = NULL, minCladeSize = 1, nSplits = 10, 
restrictNode = NULL, lowerBound = NULL, upperBound = NULL, 
tol = NULL, n.cores=1) 
{
    bounds <- matrix(c(1e-08, 1, 1e-08, 1, 1e-08, 5, 1e-08, 20, 
					   0, 1, 1e-08, 1000), 6, 2, byrow = TRUE)
    rownames(bounds) <- c("kappa", "lambda", "delta", "alpha", 
						  "psi", "rate")
	
	
	
    switch(model, bm = {
		   phy <- transformPhylo(phy = phy, model = "bm", y = y)
		   out <- likTraitPhylo(y, phy)
		   names(out) <- c("brownianVariance", "logLikelihood")
		   }, 
		   
		   kappa = {
		   kappa <- 1
		   if (is.null(lowerBound)) {
		   lowerBound <- bounds["kappa", 1]
		   }
		   if (is.null(upperBound)) {
		   upperBound <- bounds["kappa", 2]
		   }
		   
		   var.funkappa <- function(kappa) {
		   return(transformPhylo.ll(y = y, phy = phy, kappa = kappa, model = "kappa")[[2]])
		   }
		   
		   vo <- optim(kappa, var.funkappa, method = "L-BFGS-B", 
					   lower = lowerBound, upper = upperBound, control = c(fnscale = -1))
		   
		   if (modelCIs == TRUE) {
		   
		   kappa.fun <- function(param) {
		   ll <- transformPhylo.ll(y = y, phy = phy, kappa = param, model = "kappa")$logLikelihood
		   return(ll - vo$value + 1.92)
		   }
		   if (kappa.fun(lowerBound) < 0) {
		   LCI <- uniroot(kappa.fun, interval = c(lowerBound, 
											vo$par))$root
		   } else {
		   LCI <- NA
		   }
		   if (kappa.fun(upperBound) < 0) {
		   UCI <- uniroot(kappa.fun, interval = c(vo$par, upperBound))$root
		   } else {
		   UCI <- NA
		   }
		   out <- vector(mode = "list", length = 3)
		   names(out) <- c("MaximumLikelihood", "Kappa", "brownianVariance")
		   out$Kappa <- matrix(NA, 1, 3, byrow = TRUE)
		   colnames(out$Kappa) <- c("MLKappa", "LowerCI", "UpperCI")
		   out$MaximumLikelihood <- vo$value
		   out$Kappa[1, ] <- c(vo$par, LCI, UCI)
		   out$brownianVariance <- transformPhylo.ll(y=y, phy=phy, model="kappa", kappa= vo$par)$brownianVariance 
		   if (any(is.na(out$Kappa[1, 2:3]))) {
		   warning("Confidence limits fall outside the current parameter bounds - consider changing lowerBound and/or upperBound")
		   }
		   } else {
		   out <- vector(mode = "list", length = 3)
		   names(out) <- c("MaximumLikelihood", "Kappa", "brownianVariance")
		   out$Kappa <- matrix(NA, 1, 1, byrow = TRUE)
		   colnames(out$Kappa) <- c("MLKappa")
		   out$MaximumLikelihood <- vo$value
		   out$Kappa[1, ] <- c(vo$par)
		   out$brownianVariance <- transformPhylo.ll(y=y, phy=phy, model="kappa", kappa= vo$par)$brownianVariance  
		   }
		   }, 
		   
		   
		   
		   
		   lambda = {
		   lambda <- 1
		   if (is.null(lowerBound)) {
		   lowerBound <- bounds["lambda", 1]
		   }
		   if (is.null(upperBound)) {
		   upperBound <- bounds["lambda", 2]
		   }
		   
		   
		   var.funlambda <- function(lambda) {
		   
		   return(transformPhylo.ll(y = y, phy = phy, lambda = lambda, model = "lambda")[[2]])
		   }
		   
		   
		   vo <- optim(lambda, var.funlambda, method = "L-BFGS-B", 
					   lower = lowerBound, upper = upperBound, control = c(fnscale = -1))
		   
		   if (modelCIs == TRUE) {
		   lambda.fun <- function(param) {
		   ll <- transformPhylo.ll(y, phy, model = "lambda", lambda = param)$logLikelihood
		   return(ll - vo$value + 1.92)
		   }
		   if (lambda.fun(lowerBound) < 0) {
		   LCI <- uniroot(lambda.fun, interval = c(lowerBound, vo$par))$root
		   } else {
		   LCI <- NA
		   }
		   if (lambda.fun(upperBound) < 0) {
		   UCI <- uniroot(lambda.fun, interval = c(vo$par, upperBound))$root
		   } else {
		   UCI <- NA
		   }
		   out <- vector(mode = "list", length = 3)
		   names(out) <- c("MaximumLikelihood", "Lambda", "brownianVariance")
		   out$Lambda <- matrix(NA, 1, 3, byrow = TRUE)
		   colnames(out$Lambda) <- c("MLLambda", "LowerCI", 
									 "UpperCI")
		   out$MaximumLikelihood <- vo$value
		   out$Lambda[1, ] <- c(vo$par, LCI, UCI)
		   out$brownianVariance <- transformPhylo.ll(y=y, phy=phy, model="lambda", lambda= vo$par)$brownianVariance 
		   if (any(is.na(out$Lambda[1, 2:3]))) {
		   warning("Confidence limits fall outside the current parameter bounds - consider changing lowerBound and/or upperBound")
		   }
		   } else {
		   out <- vector(mode = "list", length = 3)
		   names(out) <- c("MaximumLikelihood", "Lambda")
		   out$Lambda <- matrix(NA, 1, 1, byrow = TRUE)
		   colnames(out$Lambda) <- c("MLLambda")
		   out$MaximumLikelihood <- vo$value
		   out$Lambda[1, ] <- c(vo$par)
		   out$brownianVariance <- transformPhylo.ll(y=y, phy=phy, model="lambda", lambda= vo$par)$brownianVariance 
		   }
		   }, 
		   
		   
		   
		   delta = {
		   delta <- 1
		   if (is.null(lowerBound)) {
		   lowerBound <- bounds["delta", 1]
		   }
		   if (is.null(upperBound)) {
		   upperBound <- bounds["delta", 2]
		   }
		   var.fundelta <- function(delta) {
		   return(transformPhylo.ll(y = y, phy = phy, delta = delta, model = "delta")[[2]])
		   }
		   vo <- optim(delta, var.fundelta, method = "L-BFGS-B", 
					   lower = lowerBound, upper = upperBound, control = c(fnscale = -1))
		   if (modelCIs == TRUE) {
		   delta.fun <- function(param) {
		   ll <- transformPhylo.ll(y, phy, model = "delta", delta = param)$logLikelihood
		   return(ll - vo$value + 1.92)
		   }
		   if (delta.fun(lowerBound) < 0) {
		   LCI <- uniroot(delta.fun, interval = c(lowerBound, 
											vo$par))$root
		   } else {
		   LCI <- NA
		   }
		   if (delta.fun(upperBound) < 0) {
		   UCI <- uniroot(delta.fun, interval = c(vo$par, upperBound))$root
		   } else {
		   UCI <- NA
		   }
		   out <- vector(mode = "list", length = 3)
		   names(out) <- c("MaximumLikelihood", "Delta", "brownianVariance")
		   out$Delta <- matrix(NA, 1, 3, byrow = TRUE)
		   colnames(out$Delta) <- c("MLDelta", "LowerCI", "UpperCI")
		   out$MaximumLikelihood <- vo$value
		   out$Delta[1, ] <- c(vo$par, LCI, UCI)
		   out$brownianVariance <- transformPhylo.ll(y=y, phy=phy, model="delta", delta= vo$par)$brownianVariance 
		   if (any(is.na(out$Delta[1, 2:3]))) {
		   warning("Confidence limits fall outside the current parameter bounds - consider changing lowerBound and/or upperBound")
		   }
		   } else {
		   out <- vector(mode = "list", length = 3)
		   names(out) <- c("MaximumLikelihood", "Delta", "brownianVariance")
		   out$Delta <- matrix(NA, 1, 1, byrow = TRUE)
		   colnames(out$Delta) <- c("MLDelta")
		   out$MaximumLikelihood <- vo$value
		   out$Delta[1, ] <- c(vo$par)
		   out$brownianVariance <- transformPhylo.ll(y=y, phy=phy, model="delta", delta= vo$par)$brownianVariance 
		   }
		   }, 
		   
		   
		   OU = {
		   alpha <- 0.01
		   if (is.null(lowerBound)) {
		   lowerBound <- bounds["alpha", 1]
		   }
		   if (is.null(upperBound)) {
		   upperBound <- bounds["alpha", 2]
		   }
		   var.funOU <- function(alpha) {
		   return(transformPhylo.ll(y = y, phy = phy, alpha = alpha, model = "OU")[[2]])
		   }
		   vo <- optim(alpha, var.funOU, method = "L-BFGS-B", lower = lowerBound, 
					   upper = upperBound, control = c(fnscale = -1))
		   if (modelCIs == TRUE) {
		   ou.fun <- function(param) {
		   ll <- transformPhylo.ll(y, phy, model = "OU", alpha = param)$logLikelihood
		   return(ll - vo$value + 1.92)
		   }
		   if (ou.fun(lowerBound) < 0) {
		   LCI <- uniroot(ou.fun, interval = c(lowerBound, 
											vo$par))$root
		   } else {
		   LCI <- NA
		   }
		   if (ou.fun(upperBound) < 0) {
		   UCI <- uniroot(ou.fun, interval = c(vo$par, upperBound))$root
		   } else {
		   UCI <- NA
		   }
		   out <- vector(mode = "list", length = 3)
		   names(out) <- c("MaximumLikelihood", "Alpha", "brownianVariance")
		   out$Alpha <- matrix(NA, 1, 3, byrow = TRUE)
		   colnames(out$Alpha) <- c("MLAlpha", "LowerCI", "UpperCI")
		   out$MaximumLikelihood <- vo$value
		   out$Alpha[1, ] <- c(vo$par, LCI, UCI)
		   out$brownianVariance <- transformPhylo.ll(y=y, phy=phy, model="OU", alpha= vo$par)$brownianVariance 
		   if (any(is.na(out$Alpha[1, 2:3]))) {
		   warning("Confidence limits fall outside the current parameter bounds - consider changing lowerBound and/or upperBound")
		   }
		   } else {
		   out <- vector(mode = "list", length = 3)
		   names(out) <- c("MaximumLikelihood", "Alpha", "brownianVariance")
		   out$Alpha <- matrix(NA, 1, 1, byrow = TRUE)
		   colnames(out$Alpha) <- c("MLAlpha")
		   out$MaximumLikelihood <- vo$value
		   out$Alpha[1, ] <- c(vo$par)
		   out$brownianVariance <- transformPhylo.ll(y=y, phy=phy, model="OU", alpha= vo$par)$brownianVariance 
		   }
		   }, 
		   
		   
		   psi = {
		   psi <- 1
		   if (is.null(lowerBound)) {
		   lowerBound <- bounds["psi", 1]
		   }
		   if (is.null(upperBound)) {
		   upperBound <- bounds["psi", 2]
		   }
		   var.funpsi <- function(psi) {
		   return(transformPhylo.ll(y = y, phy = phy, psi = psi, model = "psi")[[2]])
		   }
		   vo <- optim(psi, var.funpsi, method = "L-BFGS-B", lower = lowerBound, 
					   upper = upperBound, control = c(fnscale = -1))
		   if (modelCIs == TRUE) {
		   psi.fun <- function(param) {
		   ll <- transformPhylo.ll(y, phy, model = "psi", psi = param)$logLikelihood
		   return(ll - vo$value + 1.92)
		   }
		   if (psi.fun(lowerBound) < 0) {
		   LCI <- uniroot(psi.fun, interval = c(lowerBound, 
											vo$par))$root
		   } else {
		   LCI <- NA
		   }
		   if (psi.fun(upperBound) < 0) {
		   UCI <- uniroot(psi.fun, interval = c(vo$par, upperBound))$root
		   } else {
		   UCI <- NA
		   }
		   out <- vector(mode = "list", length = 3)
		   names(out) <- c("MaximumLikelihood", "psi", "brownianVariance")
		   out$psi <- matrix(NA, 1, 3, byrow = TRUE)
		   colnames(out$psi) <- c("MLpsi", "LowerCI", "UpperCI")
		   out$MaximumLikelihood <- vo$value
		   out$psi[1, ] <- c(vo$par, LCI, UCI)
		   out$brownianVariance <- transformPhylo.ll(y=y, phy=phy, model="psi", psi= vo$par)$brownianVariance 
		   if (any(is.na(out$psi[1, 2:3]))) {
		   warning("Confidence limits fall outside the current parameter bounds - consider changing lowerBound and/or upperBound")
		   }
		   } else {
		   out <- vector(mode = "list", length = 3)
		   names(out) <- c("MaximumLikelihood", "psi", "brownianVariance")
		   out$psi <- matrix(NA, 1, 1, byrow = TRUE)
		   colnames(out$psi) <- c("MLpsi")
		   out$MaximumLikelihood <- vo$value
		   out$psi[1, ] <- c(vo$par)
		   out$brownianVariance <- transformPhylo.ll(y=y, phy=phy, model="psi", psi= vo$par)$brownianVariance 
		   }
		   }, 
		   
		   
		   free = {
		   if (is.null(lowerBound)) {
		   lowerBound <- bounds["rate", 1]
		   }
		   if (is.null(upperBound)) {
		   upperBound <- bounds["rate", 2]
		   }
		   branchRates <- rep(1, length(phy$edge.length))
		   var.funfree <- function(branchRates) {
		   return(transformPhylo.ll(y, phy, branchRates = branchRates, model = "free")[[2]])
		   }
		   vo <- optim(branchRates, var.funfree, method = "L-BFGS-B", 
					   lower = 0, control = c(fnscale = -1, maxit = 10, 
											  factr = 1e+14))
		   phy2 <- phy
		   phy2$edge.length <- phy$edge.length * vo$par
		   out <- vector(mode = "list", length = 3)
		   names(out) <- c("MaximumLikelihood", "Rates", "Convergence")
		   out$MaximumLikelihood <- transformPhylo.ML(y, phy = phy2, model = "bm")[[2]]
		   out$Rates <- vo$par
		   if (vo$convergence == 0) {
		   out$Convergence <- "Successful"
		   } else {
		   out$Convergence <- "Failed"
		   }
		   }, 
		   
		   
		   clade = {
		   if (is.null(lowerBound)) {
		   lowerBound <- bounds["rate", 1]
		   }
		   if (is.null(upperBound)) {
		   upperBound <- bounds["rate", 2]
		   }
		   cladeRates <- rep(1, length(nodeIDs))
		   var.funclade <- function(cladeRates) {
		   return(transformPhylo.ll(y, phy, nodeIDs = nodeIDs, 
									cladeRates = cladeRates, model = "clade", rateType = rateType)[[2]])
		   }
		   vo <- optim(cladeRates, var.funclade, method = "L-BFGS-B", 
					   lower = lowerBound, upper = upperBound, control = c(fnscale = -1))
		   if (modelCIs == TRUE) {
		   free.fun <- function(param) {
		   ll <- transformPhylo.ll(y, phyClade, model = "clade", 
								   nodeIDs = SingleNode, rateType = whichRateType, cladeRates = param)$logLikelihood
		   return(ll - vo$value + 1.92)
		   }
		   phyClade <- transformPhylo(phy, model = "clade", 
									  nodeIDs = nodeIDs, rateType = rateType, cladeRates = vo$par, y = y)
		   out <- vector(mode = "list", length = 2)
		   names(out) <- c("MaximumLikelihood", "Rates")
		   out$Rates <- matrix(NA, length(nodeIDs), 4, byrow = TRUE)
		   colnames(out$Rates) <- c("node", "MLRate", "LowerCI", 
									"UpperCI")
		   out$MaximumLikelihood <- vo$value
		   for (i in 1:length(nodeIDs)) {
		   SingleNode <- nodeIDs[i]
		   whichRateType <- rateType[i]
		   phyClade <- transformPhylo(phyClade, model = "clade", 
									  nodeIDs = SingleNode, rateType = whichRateType, 
									  cladeRates = 1/vo$par[i], y = y)
		   if (free.fun(lowerBound) < 0) {
		   LCI <- uniroot(free.fun, interval = c(lowerBound, 
											vo$par[i]))$root
		   } else {
		   LCI <- NA
		   }
		   if (free.fun(upperBound) < 0) {
		   UCI <- uniroot(free.fun, interval = c(vo$par[i], 
											upperBound))$root
		   } else {
		   UCI <- NA
		   }
		   out$Rates[i, ] <- c(nodeIDs[i], vo$par[i], LCI, 
							   UCI)
		   }
		   if (any(is.na(out$Rates[, 3:4]))) {
		   warning("Confidence limits fall outside the current parameter bounds - consider changing lowerBound and/or upperBound")
		   }
		   } else {
            phyClade <- transformPhylo(phy, model = "clade", nodeIDs = nodeIDs, rateType = rateType, cladeRates = vo$par, y = y)
            out <- vector(mode = "list", length = 3)
		   names(out) <- c("MaximumLikelihood", "Rates", "brownianVariance")
		   out$Rates <- matrix(NA, length(nodeIDs), 2, byrow = TRUE)
		   colnames(out$Rates) <- c("node", "MLRate")
		   out$MaximumLikelihood <- vo$value
		   out$Rates[, 1] <- nodeIDs
		   out$Rates[, 2] <- vo$par
		   out$brownianVariance <- transformPhylo.ll(y=y, phy=phyClade, nodeIDs = nodeIDs, rateType = rateType, model="clade", cladeRates= vo$par)$brownianVariance
		   }
		   }, 
		   
		   
		   tm1 = {  
		   if (is.null(lowerBound)) {
		   lowerBound <- bounds["rate", 1]
		   }
		   if (is.null(upperBound)) {
		   upperBound <- bounds["rate", 2]
		   }
		   nodes <- unique(sort(phy$edge[, 2]))
		   nodeDepth <- node.depth(phy)
		   nodesByRichness <- cbind(richness = nodeDepth[nodes], 
									node = nodes)
		   searchNode <- nodesByRichness[nodesByRichness[, 1] >= 
		   minCladeSize, 2]
		   if (is.null(restrictNode) == FALSE) {
		   ngroups <- length(restrictNode)
		   cm <- clade.matrix(phy)
		   colnames(cm$clade.matrix) <- cm$tip.label
		   cmMat <- cm$clade.matrix
		   skipNodes <- numeric()
		   for (i in 1:ngroups) {
		   matDrop <- cmMat[, restrictNode[[i]]]
		   nodeDrop <- rowSums(matDrop) < length(restrictNode[[i]]) & 
		   rowSums(matDrop) >= 1
		   skipNodes <- c(skipNodes, as.numeric(rownames(matDrop)[nodeDrop]))
		   }
		   skipNodes <- unique(c(skipNodes, phy$edge[phy$edge.length < 
								 tol, 2]))
		   searchNode <- setdiff(searchNode, skipNodes)
		   } else {
		   skipNodes <- unique(phy$edge[phy$edge.length < tol, 
							   2])
		   searchNode <- setdiff(searchNode, skipNodes)
		   }
		   n <- length(y)
		   BERateOut <- matrix(NA, nrow = nSplits, ncol = (6 + nSplits))
		   
		   
		   fullModelOut <- vector(mode = "list", length = nSplits)
		   
		   for (k in 1:nSplits) {
		   cladeMembers <- matrix(NA, ncol = k, nrow = length(phy$edge[, 
															  1]))
		   if (k == 1) {
		   searchNode <- searchNode
		   } else {
		   searchNode <- searchNode[!(searchNode %in% bestNodes)]
		   }
		   for (i in searchNode) {
		   if (k == 1) {
		   currentNodeIDs <- i
		   } else {
		   currentNodeIDs <- c(bestNodes, i)
		   }
		   cladeRates <- rep(1, length(currentNodeIDs))
		   var.funclade <- function(cladeRates) {
		   return(transformPhylo.ll(y, phy, nodeIDs = currentNodeIDs, 
									rateType = rep("clade", length(currentNodeIDs)), 
									cladeRates = cladeRates, model = "clade")[[2]])
		   }
		   currentCladeModel <- optim(cladeRates, var.funclade, 
									  method = "L-BFGS-B", lower = lowerBound, upper = upperBound, 
									  control = c(fnscale = -1))
		   fullModelOut[[k]] <- rbind(fullModelOut[[k]], 
									  c(node = as.integer(i), shiftPos = 1, ML = currentCladeModel$value, 
										currentCladeModel$par))
		   currentModel <- currentCladeModel
		   shiftPos = 1
		   param <- (k * 2) + 1
		   currentModel <- list(currentModel, i, cladeMembers)
		   currentML <- currentModel[[1]]$value
		   AIC <- -2 * currentModel[[1]]$value + 2 * param
		   AICc <- -2 * currentModel[[1]]$value + 2 * param + 
		   ((2 * param * (param + 1))/(n - param - 1))
		   if (i == min(searchNode)) {
		   BERateOut[k, 1:(6 + k)] <- c(currentModel[[2]], 
										shiftPos, currentModel[[1]]$value, param, 
										AIC, AICc, currentModel[[1]]$par)
		   } else {
		   if (currentML > BERateOut[k, 3]) {
		   BERateOut[k, 1:(6 + k)] <- c(currentModel[[2]], 
										shiftPos, currentModel[[1]]$value, param, AIC, AICc, currentModel[[1]]$par)
		   }
		   }
		   }
		   if (k == 1) {
		   bestNodes <- BERateOut[k, 1]
		   } else {
		   bestNodes <- c(bestNodes, BERateOut[k, 1])
		   }
		   print(BERateOut)
		   }
		   BERateSummary <- matrix(NA, ncol = (6 + nSplits), nrow = 1)
		   BERateSummary <- as.data.frame(BERateSummary)
		   MLsingle <- likTraitPhylo(y, phy)[[2]]
		   AICsingle <- -2 * MLsingle + 2 * 2
		   AICcsingle <- -2 * MLsingle + 2 * 2 + ((2 * 2 * (2 + 
															1))/(n - 2 - 1))
		   BERateSummary[1, 1:6] <- c(0, 1, MLsingle, 2, AICsingle, 
									  AICcsingle)
		   BERateSummary <- rbind(BERateSummary, BERateOut)
		   colnames(BERateSummary) <- c("node", "shiftPos", "ML", 
										"k", "AIC", "AICc", c(BERateSummary[2:nrow(BERateSummary), 
															  1]))
		   if (sum(BERateSummary[, "shiftPos"] == 1) > 0) {
		   BERateSummary[which(BERateSummary[, "shiftPos"] == 
							   1), "shiftPos"] <- "clade"
		   }
		   out <- list(as.data.frame(BERateSummary), fullModelOut, 
					   y, phy)
		   }, 
		   
		   
		   
		   
		   
		   
		   
		   
		   
		   tm2 = {
		   if (is.null(lowerBound)) {
		   lowerBound <- bounds["rate", 1]
		   }
		   if (is.null(upperBound)) {
		   upperBound <- bounds["rate", 2]
		   }
		   nodes <- unique(sort(phy$edge[, 2]))
		   nodeDepth <- node.depth(phy)
		   nodesByRichness <- cbind(richness = nodeDepth[nodes], 
									node = nodes)
		   searchNode <- nodesByRichness[nodesByRichness[, 1] >= 
		   minCladeSize, 2]
		   searchNodeType <- c(searchNode, searchNode[searchNode > 
							   Ntip(phy)])
		   searchType <- c(rep("clade", length(searchNode)), rep("branch", 
																 sum(searchNode > Ntip(phy))))
		   searchNodeTypeAll <- data.frame(searchType, searchNodeType, 
										   stringsAsFactors = FALSE)
		   if (is.null(restrictNode) == FALSE) {
		   ngroups <- length(restrictNode)
		   cm <- clade.matrix(phy)
		   colnames(cm$clade.matrix) <- cm$tip.label
		   cmMat <- cm$clade.matrix
		   skipNodes <- numeric()
		   for (i in 1:ngroups) {
		   matDrop <- cmMat[, restrictNode[[i]]]
		   nodeDrop <- rowSums(matDrop) < length(restrictNode[[i]]) & 
		   rowSums(matDrop) >= 1
		   skipNodes <- c(skipNodes, as.numeric(rownames(matDrop)[nodeDrop]))
		   }
		   skipNodes <- unique(c(skipNodes, phy$edge[phy$edge.length < 
								 tol, 2]))
		   searchNode <- setdiff(searchNode, skipNodes)
		   } else {
		   skipNodes <- unique(phy$edge[phy$edge.length < tol, 
							   2])
		   searchNode <- setdiff(searchNode, skipNodes)
		   }
		   tm2.fun <- function(whichSearchNode) {
		   which(searchNodeTypeAll[, 2] == whichSearchNode)
		   }
		   searchNodeTypeAll <- searchNodeTypeAll[sort(unlist(sapply(searchNode, 
																	 tm2.fun))), ]
		   n <- length(y)
		   BERateOut <- matrix(NA, nrow = nSplits, ncol = (6 + (nSplits)))
		   fullModelOut <- vector(mode = "list", length = nSplits)
		   bestModel <- matrix(NA, ncol = 2, nrow = nSplits)
		   
		   cladeMembersObj <- allCladeMembers(phy)
		   
		   fullModelOut <- vector(mode = "list", length = nSplits)
		   
		   for (k in 1:nSplits) {
		   cladeMembers <- matrix(NA, ncol = k, nrow = length(phy$edge[, 
															  1]))
		   bestModel[bestModel[, 2] == "1", 2] <- "clade"
		   bestModel[bestModel[, 2] == "2", 2] <- "branch"
		   if (k == 1) {
		   searchNodeTypeAll <- searchNodeTypeAll
		   } else {
		   dropNodeRow <- which(searchNodeTypeAll[, 1] == 
								bestModel[k - 1, 2] & searchNodeTypeAll[, 2] == 
								as.numeric(bestModel[k - 1, 1]))
		   rowDropLogical <- row(searchNodeTypeAll)[, 1] != dropNodeRow
		   searchNodeTypeAll <- searchNodeTypeAll[rowDropLogical, ]
		   }
           
### Weird stuff going on here. Updating of currentNodeIDs seems to be the problem.  Needs to be taken out of the function? Compare old with new version to see where problem is.
		   
		   if (k == 1) {
		   currentNodeIDs <- NULL
		   } else {
		   currentNodeIDs <- na.omit(as.numeric(bestModel[, 1]))
		   }
		   if (k == 1) {
		   currentRateType <- NULL
		   } else {
		   currentRateType <- na.omit(bestModel[, 2])
		   }
		   
		   
		   foo2 <- function(x, k, currentNodeIDs=currentNodeIDs, currentRateType=currentRateType, cladeMembersObj=cladeMembersObj){
		   
		   currentNodeIDs <- c(currentNodeIDs, x[,2])
		   currentRateType <- c(currentRateType, x[,1])
		   
		   
		   cladeRates <- rep(1, length(currentNodeIDs))
		   var.funclade <- function(cladeRates) {
		   return(transformPhylo.ll(y, phy, nodeIDs = currentNodeIDs, 
									rateType = currentRateType, cladeRates = cladeRates, cladeMembersObj=cladeMembersObj,
									model = "clade")[[2]])
		   }
		   currentCladeModel <- optim(cladeRates, var.funclade, 
									  method = "L-BFGS-B", lower = lowerBound, upper = upperBound, 
									  control = c(fnscale = -1))
		   
		   if (x[,1] == "clade") {
		   shiftPos <- 1
		   }
		   if (x[,1] == "branch") {
		   shiftPos <- 2
		   }

		   currentModel <- currentCladeModel
		   param <- k*2 + 1
		   currentModel <- list(currentModel, x[2], cladeMembers)
		   currentML <- currentModel[[1]]$value
		   AIC <- -2 * currentModel[[1]]$value + 2 * param
		   AICc <- -2 * currentModel[[1]]$value + 2 * param + 
		   ((2 * param * (param + 1))/(n - param - 1))
		   
		   
		   out <- c(currentModel[[2]][,1], 
                    shiftPos, currentModel[[1]]$value, param, 
                    AIC, AICc, currentModel[[1]]$par)
		   
		   return(out)}
		   		   
		   searchNodeList <- vector(mode="list", length=length(searchNodeTypeAll[,1]))
		   for (i in 1:length(searchNodeTypeAll[,1])) {searchNodeList[[i]] <- searchNodeTypeAll[i,]}
		   		   
		   all.output <- mclapply(searchNodeList, FUN=foo2, k=k, currentNodeIDs=currentNodeIDs, currentRateType=currentRateType, cladeMembersObj=cladeMembersObj, mc.cores=n.cores)
		   all.output <- unlist(all.output)
		   all.output <- matrix(all.output, nrow=length(searchNodeList) , byrow=TRUE)
		   best.out <- all.output[which(all.output[,3]==max(all.output[,3])),]
		   
		   fullModelOut[[k]] <- all.output
		   
			if ((length(BERateOut[k,])-length(best.out)) > 0) {
		   BERateOut[k,] <- c(best.out, rep(NA, length(BERateOut[k,])-length(best.out)))
		   } else {BERateOut[k,] <- best.out}
		   
		   
		   if (k == 1) {
		   bestModel[1, ] <- BERateOut[k, 1:2]
		   } else {
		   bestModel[k, ] <- BERateOut[k, 1:2]
		   }
		   print(BERateOut)
		   }
		   
		   BERateSummary <- matrix(NA, ncol = (6 + nSplits), nrow = 1)
		   BERateSummary <- as.data.frame(BERateSummary)
		   MLsingle <- likTraitPhylo(y, phy)[[2]]
		   AICsingle <- -2 * MLsingle + 2 * 2
		   AICcsingle <- -2 * MLsingle + 2 * 2 + ((2 * 2 * (2 + 
															1))/(n - 2 - 1))
		   BERateSummary[1, 1:6] <- c(0, 1, MLsingle, 2, AICsingle, 
									  AICcsingle)
		   BERateSummary <- rbind(BERateSummary, BERateOut)
		   colnames(BERateSummary) <- c("node", "shiftPos", "ML", 
										"k", "AIC", "AICc", c(BERateSummary[2:nrow(BERateSummary), 
															  1]))
		   if (sum(BERateSummary[, "shiftPos"] == 1) > 0) {
		   BERateSummary[which(BERateSummary[, "shiftPos"] == 
							   1), "shiftPos"] <- "clade"
		   }
		   if (sum(BERateSummary[, "shiftPos"] == 2) > 0) {
		   BERateSummary[which(BERateSummary[, "shiftPos"] == 
							   2), "shiftPos"] <- "branch"
		   }
		   out <- list(as.data.frame(BERateSummary), fullModelOut, 
					   y, phy)
		   })
    return(out)
}
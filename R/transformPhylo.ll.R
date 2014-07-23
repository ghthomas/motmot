transformPhylo.ll<-function (y = NULL, phy, model = NULL, meserr = NULL, 	kappa = NULL, lambda = NULL, delta = NULL, alpha = NULL, psi = NULL, la = NULL, nodeIDs = NULL, rateType = NULL, branchRates = NULL, cladeRates = NULL, branchLabels = NULL, twosigma = FALSE, covPIC = TRUE)
{
    switch(model, bm = {
        transformPhy <- transformPhylo(phy = phy, model = "bm",
        meserr = meserr, y = y)
    }, kappa = {
        transformPhy <- transformPhylo(phy = phy, model = "kappa",
        kappa = kappa, meserr = meserr, y = y)
    }, lambda = {
        transformPhy <- transformPhylo(phy = phy, model = "lambda",
        lambda = lambda, meserr = meserr, y = y)
    }, delta = {
        transformPhy <- transformPhylo(phy = phy, model = "delta",
        delta = delta, meserr = meserr, y = y)
    }, free = {
        transformPhy <- transformPhylo(phy = phy, model = "free",
        branchRates = branchRates, meserr = meserr, y = y)
    }, clade = {
        transformPhy <- transformPhylo(phy = phy, model = "clade",
        nodeIDs = nodeIDs, cladeRates = cladeRates, rateType = rateType,
        meserr = meserr, y = y)
    }, OU = {
        transformPhy <- transformPhylo(phy = phy, model = "OU",
        alpha = alpha, meserr = meserr, y = y)
    }, psi = {
        transformPhy <- transformPhylo(phy = phy, model = "psi",
        psi = psi, meserr = meserr, y = y, la = la)
    }, multipsi = {
        transformPhy <- transformPhylo(phy = phy, branchLabels = branchLabels, model = "multipsi", psi = psi, meserr = meserr, y = y, la = la)
    })
    return(likTraitPhylo(y = y, phy = transformPhy, covPIC = covPIC))
}

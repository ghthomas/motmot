transformPhylo.ll <- function(y, phy, model=NULL, meserr=NULL, kappa=NULL, lambda=NULL, delta=NULL, alpha=NULL, psi=NULL, nodeIDs=NULL, rateType=NULL, branchRates=NULL, cladeRates=NULL, cladeMembersObj=NULL) {
	
	switch(model,		  
		   
		   "bm" = {
		   transformPhy <- transformPhylo(phy=phy, model="bm", y=y)
		   },
		   
		   "kappa" = {
		   transformPhy <- transformPhylo(phy=phy, model="kappa", kappa=kappa, y=y)
		   },
		   
		   "lambda" = {
		   transformPhy <- transformPhylo(phy=phy, model="lambda", lambda=lambda, y=y)
		   },
		   
		   "delta" = {
		   transformPhy <- transformPhylo(phy=phy, model="delta", delta=delta, y=y)
		   },
		   
		   "free" = {
		   transformPhy <- transformPhylo(phy=phy, model="free", branchRates=branchRates, y=y)
		   },
		   
		   "clade" = {
		   transformPhy <- transformPhylo(phy=phy, model="clade", nodeIDs=nodeIDs, cladeRates=cladeRates, rateType=rateType, y=y, cladeMembersObj=cladeMembersObj)
		   },
		   
		   "OU" = {
		   transformPhy <- transformPhylo(phy=phy, model="OU", alpha=alpha, y=y)
		   },
		   
		   "psi" = {
		   transformPhy <- transformPhylo(phy=phy, model="psi", psi=psi, y=y)
		   }
		   
		   )
	
	return(likTraitPhylo(y=y, phy=transformPhy))
}

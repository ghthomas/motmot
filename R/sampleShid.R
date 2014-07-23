sampleShid <- function(phy, la = NULL, mu = NULL, useMean = FALSE)
{
	if(is.null(la)|is.null(mu))
    stop("Please provide values for lambda and mu ")
	phy$node.label<-NULL
	br<-branching.times(phy)
	to<-br[as.character(phy$edge[,1])]
	start<-br[as.character(phy$edge[,1])]
	te<-start-phy$edge.length
	te[te<0]<-0
	if(mu > 0){
		if(la == mu){
			expSh <- 2 * la * (to - te) + 2 * log((1 + la * te) / (1 + la * to))
		}else{
			expSh <- 2 * la * (to - te) + 2 * log((la * exp(te * (la - mu)) - mu) / (la * exp(to * (la - mu)) - mu))
		}
		if(useMean){
			Shid <- expSh
		}else{
			Shid <- rpois(length(expSh), expSh)
		}
	}else{
		Shid <- rep(0, length(to))
	}
	phy$Shid <- Shid
	phy
}
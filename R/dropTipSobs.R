dropTipSobs<-function(phy, tip)
{
	phy1 <- phy
	phy1$edge.length <- phy$edge.length^0
	phy1 <- drop.tip(phy1, tip)
	phy2 <- drop.tip(phy, tip)
	phy2$Sobs <- phy1$edge.length
	phy2
}

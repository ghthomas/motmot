\name{transformPhylo.ll}
\alias{transformPhylo.ll}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log-likelhood for models of trait evoluion. 
}
\description{
Fits likelihood models for various models of continuous character evolution.}
\usage{
transformPhylo.ll(y, phy, model = NULL, meserr=NULL, kappa = NULL, lambda = NULL, delta = NULL, alpha = NULL, psi=NULL, nodeIDs = NULL, rateType = NULL, branchRates = NULL, cladeRates = NULL, cladeMembersObj=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{y}{
A matrix of trait values.
}
  \item{phy}{
An object of class "phylo" (see ape package).
}
 \item{model}{
The model of trait evolution (see details).}
  \item{meserr}{
A vector (or matrix) of measurement error for each tip. This is only applicable to univariate analyses.}
  \item{kappa}{
Value of kappa transform.}
  \item{lambda}{
Value of lambda transform.}
  \item{delta}{
Value of delta transform.}
  \item{alpha}{
Value of alpha (OU) transform.}
  \item{psi}{
Value of psi transform.}
  \item{nodeIDs}{
Integer - ancestral nodes of clades.}
  \item{rateType}{
  If model="clade", a vector specifying if rate shift occurs in a clade ("clade") or on the single branch leading to a clade ("branch").}
\item{branchRates}{
Numeric vector specifying relative rates for individual branches.}
  \item{cladeRates}{
Numeric vector specifying telative rates for clades.}
  \item{cladeMembersObj}{
Matrix of clade membership}
}
\details{
This function fits likelihood models (see below) for continuous character evolution where the parameter values are set a priori. The function returns the log-likihood and the Brownian variance (or variance covariance matrix).

model="bm"- Brownian motion (constant rates random walk)

model="kappa" - fits Pagel's kappa by raising all branch lengths to the power kappa. As kappa approaches zero, trait change becomes focused at branching events. For complete phylogenies, if kappa approaches zero this infers speciational trait change. 

model="lambda" - fits Pagel's lambda to estimate phylogenetic signal by multiplying all internal branches of the tree by lambda, leaving tip branches as their original length (root to tip distances are unchanged);

model="delta" - fits Pagel's delta by raising all node depths to the power delta. If delta <1, trait evolution is concentrated early in the tree whereas if delta >1 trait evolution is concentrated towards the tips. Values of delta above one can be difficult to fit reliably.

model="free" - fits Mooer's et al's free model where each branch has its own rate of trait evolution. This can be a useful exploratory analysis but it is slow due to the number of parameters, particularly for large trees.

model="clade" - fits a model where particular clades are a priori hypothesised to have different rates of trait evolution (see O'Meara et al. 2006; Thomas et al. 2006, 2009). Clades are specified using nodeIDs and are defined as the mrca node. Unique rates for each clade are specified using cladeRates. rateType specifies whether the rate shift occurs in the stem clade or on the single branch leading to the clade.

model="OU" - fits an Ornstein-Uhlenbeck model - a random walk with a central tendency proportional to alpha. High values of alpha can be interpreted as evidence of evolutionary constraints, stabilising selection or weak phylogenetic signal.

model="psi" - fits a model to assess to the relative contributions of speciation and gradual evolution to a trait's evolutionary rate (Ingram 2010). 
}

\value{
  \item{brownianVariance }{Brownian variance (or covariance for multiple traits) given the data and phylogeny}
  \item{logLikelihood }{The log-likelihood of the model and data}
  


}
\references{

Felsenstein J. 1973. Maximum-likelihood estimation of evolutionary trees from continuous characters. Am. J. Hum. Genet. 25, 471-492.

Felsenstein J. 1985. Phylogenies and the comparative method. American Naturalist 125, 1-15.

Freckleton RP & Jetz W. 2009. Space versus phylogeny: disentangling phylogenetic and spatial signals in comparative data. Proc. Roy. Soc. B 276, 21-30. 

Ingram T. 2010. Speciation along a depth gradient in a marine adaptive radiation. Proceeding of the Royal Society B. In press.

Mooers AO, Vamosi S, & Schluter D. 1999. Using phylogenies to test macroevolutionary models of trait evolution: sexual selection and speciation in Cranes (Gruinae). American Naturalist 154, 249-259.

O'Meara BC, Ane C, Sanderson MJ & Wainwright PC. 2006. Testing for different rates of continuous trait evolution using likelihood. Evolution 60, 922-933

Pagel M. 1997. Inferring evolutionary processes from phylogenies. Zoologica Scripta 26, 331-348.

Pagel M. 1999 Inferring the historical patterns of biological evolution. Nature 401, 877-884.

Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments and island effects. Evolution 63, 2017-2030.
}
\author{
Gavin Thomas}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
# Data and phylogeny
data(anolis.tree)
data(anolis.data)

# anolis.data is not matrix and contains missing data so put together matrix of # relevant traits (here female and male snout vent lengths) and remove species 
# with missing data from the matrix and phylogeny
anolisSVL <- data.matrix(anolis.data)[,c(5,6)]
anolisSVL[,1] <- log(anolisSVL[,1])
anolisSVL[,2] <- log(anolisSVL[,2])

tree <- drop.tip(anolis.tree, names(attr(na.omit(anolisSVL), "na.action")))
anolisSVL <- na.omit(anolisSVL)

# log likelihood of kappa = 0.1 or 1
transformPhylo.ll(anolisSVL, phy=tree, model="kappa", kappa=0.1)
transformPhylo.ll(anolisSVL, phy=tree, model="kappa", kappa=1)

# log likelihood of lambda = 0.01 or 1
transformPhylo.ll(anolisSVL, phy=tree, model="lambda", lambda=0.01)
transformPhylo.ll(anolisSVL, phy=tree, model="lambda", lambda=1)

# log likelihood of delta = 1.5 or 1
transformPhylo.ll(anolisSVL, phy=tree, model="delta", delta=1.5)
transformPhylo.ll(anolisSVL, phy=tree, model="delta", delta=1)

# log likelihood of alpha = 0.001 or 2
transformPhylo.ll(anolisSVL, phy=tree, model="OU", alpha=0.001)
transformPhylo.ll(anolisSVL, phy=tree, model="OU", alpha=2)

# log likelihood of psi = 0 (gradual) or 1 (speciational)
transformPhylo.ll(anolisSVL, phy=tree, model="psi", psi=0)
transformPhylo.ll(anolisSVL, phy=tree, model="psi", psi=1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

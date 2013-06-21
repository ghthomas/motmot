\name{transformPhylo}
\alias{transformPhylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Phylogenetic tree transformations
}
\description{
Transforms the branch lengths of a phylo object according to a model of trait evolution (see details).}
\usage{
transformPhylo(phy, model = NULL, y=NULL, kappa = NULL, lambda = NULL, delta = NULL, alpha=NULL, psi=NULL, nodeIDs = NULL, rateType=NULL, branchRates=NULL, cladeRates = NULL, cladeMembersObj=NULL, meserr=NULL, sigmaScale=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{phy}{
An object of class "phylo" (see ape package).
}
 \item{model}{
The model of trait evolution (see details).}
\item{y}{
A matrix of trait values.}
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
  \item{meserr}{
A vector (or matrix) of measurement error for each tip. This is only applicable to univariate analyses.}
  \item{sigmaScale}{
Scalar for measurement error.}
}

\details{
Transforms the branch lengths of a phylo object according to one of the following models:

model="bm"- Brownian motion (constant rates random walk)

model="kappa" - fits Pagel's kappa by raising all branch lengths to the power kappa. As kappa approaches zero, trait change becomes focused at branching events. For complete phylogenies, if kappa approaches zero this infers speciational trait change. 

model="lambda" - fits Pagel's lambda to estimate phylogenetic signal by multiplying all internal branches of the tree by lambda, leaving tip branches as their original length (root to tip distances are unchanged);

model="delta" - fits Pagel's delta by raising all node depths to the power delta. If delta <1, trait evolution is concentrated early in the tree whereas if delta >1 trait evolution is concentrated towards the tips. Values of delta above one can be difficult to fit reliably.

model="free" - fits Mooer's et al's (1999) free model where each branch has its own rate of trait evolution. This can be a useful exploratory analysis but it is slow due to the number of parameters, particularly for large trees.

model="clade" - fits a model where particular clades are a priori hypothesised to have different rates of trait evolution (see O'Meara et al. 2006; Thomas et al. 2006, 2009). Clades are specified using nodeIDs and are defined as the mrca node. Unique rates for each clade are specified using cladeRates. rateType specifies whether the rate shift occurs in the stem clade or on the single branch leading to the clade.

model="OU" - fits an Ornstein-Uhlenbeck model - a random walk with a central tendency proportional to alpha. High values of alpha can be interpreted as evidence of evolutionary constraints, stabilising selection or weak phylogenetic signal.

model="psi" - fits a model to assess to the relative contributions of speciation and gradual evolution to a trait's evolutionary rate (Ingram 2010). 
}

\value{
  \item{phy }{A phylo object}

}
\references{

Ingram T. 2010. Speciation along a depth gradient in a marine adaptive radiation. Proceeding of the Royal Society B. In press.

Mooers AO, Vamosi S, & Schluter D. 1999. Using phylogenies to test macroevolutionary models of trait evolution: sexual selection and speciation in Cranes (Gruinae). American Naturalist 154, 249-259.

O'Meara BC, Ane C, Sanderson MJ & Wainwright PC. 2006. Testing for different rates of continuous trait evolution using likelihood. Evolution 60, 922-933

Pagel M. 1997. Inferring evolutionary processes from phylogenies. Zoologica Scripta 26, 331-348.

Pagel M. 1999 Inferring the historical patterns of biological evolution. Nature 401, 877-884.

Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments and island effects. Evolution 63, 2017-2030.}
\author{
Gavin Thomas}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

# Apply delta transform to anolis tree
data(anolis.tree)
anolis.treeDelta <- transformPhylo(phy=anolis.tree, model="delta", delta=0.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{transformPhylo.sim}
\alias{transformPhylo.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate trait data
}
\description{
Simulates continuous trait data under a range of evolutionary models}
\usage{
transformPhylo.sim(phy, n = 1, x = NULL, model = NULL, kappa = NULL, lambda = NULL, delta = NULL, alpha = NULL, psi = NULL, nodeIDs = NULL, rateType = NULL, cladeRates = NULL, branchRates=NULL, rate = NULL, group.means = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
An object of class "phylo" (see ape package).
}
  \item{n}{
Number of simulations
}
  \item{x}{
Vector, matrix or data.frame (with taxon names as names or rownames) of categories for each species. Only applicable if model="mixedRate" 
}
  \item{model}{
The model of trait evolution (see details).
}
  \item{kappa}{
Value of kappa transform.
}
  \item{lambda}{
Value of lambda transform.
}
  \item{delta}{
Value of delta transform.
}
  \item{alpha}{
Value of alpha (OU) transform.
}
  \item{psi}{
Value of psi transform.
}
  \item{nodeIDs}{
Integer - ancestral nodes of clades.
}
  \item{rateType}{
If rate is to be inferred on single branches leading to clades then rateType="branches". If rate is inferred for whole clade, rateType="clade"}
  \item{branchRates}{
Numeric vector specifying relative rates for individual branches.}
\item{cladeRates}{
Numeric vector specifying relative rates for clades.
}
  \item{rate}{
a vector of relative rate parameters. The length of the vector is equal to the number of rates being estimated. 
}
  \item{group.means}{
a vector of the relative difference in means between rate categories, expressed as a scalar applied 
to the expected standard deviation (see Ricklefs 2006)}
}

\value{
Returns a matrix of simulated dated with taxon names as rownames (number of columns=n).
}

\references{
Ricklefs RE. 2006. Time, species, and the generation of trait variation in clades. Systematic 
Biology 55, 151-159.

Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments
and island effects. Evolution 63, 2017-2030.}

\author{
Gavin Thomas}

\examples{
data(anolis.tree)
data(anolis.data)

# Simulate 10 sets of data with kappa=0.1 using the anolis tree
sim.dat1 <- transformPhylo.sim(phy=anolis.tree, n=10, model="kappa", kappa=0.1)

# Simulate 10 sets of data where rates and means differ between to the categories defined by "x"
x <- anolis.data$geo_ecomorph
names(x) <-  rownames(anolis.data)
sim.dat2 <- transformPhylo.sim(phy=anolis.tree, n=10, x=x, model="mixedRate", rate=c(1,1,2,4), group.means=c(0,5,0,0))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

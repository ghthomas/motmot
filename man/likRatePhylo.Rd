\name{likRatePhylo}
\alias{likRatePhylo}
\title{
Log-likelihood rate estimation for traits and phylogenies
}

\description{
This function calculates the log-likelihood, phylogenetic mean, and Brownian variance for a trait 
and a phylogeny transformed according to variation in relative rates. }

\usage{
likRatePhylo(rateData, rate = NULL, common.mean = FALSE, lambda.est, lambda, meserr, sigmaScale)
}

%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{rateData}{
an object of class "rateData"}

  \item{rate}{
a vector of relative rate parameters. The length of the vector is equal to the number of rates being
estimated. If \code{rate=NULL} then rates are equal.}

  \item{common.mean}{a logical specififying whether each rate category should have its own mean (\code{common.mean=FALSE})
or all categories should have the same mean (\code{common.mean=FALSE}). See Thomas et al. (2009) 
for a discussion on the impact of assumptions about mean on rate estimates.}

  \item{lambda.est}{
Logical. Fit Pagel's lambda.}

  \item{lambda}{
Numeric value for lambda from 0-1. }
  \item{meserr}{
Logical. Include measurement error. }
  \item{sigmaScale}{
Scalar for measurement error relative to tree. }
}
\value{
\item{ll }{log-likelihood of the model}
\item{mu }{phylogenetically corrected mean(s)}
\item{s2 }{Brownian variance}
}
\references{
Thomas GH, Freckleton RP, & Szekely T. 2006. Comparative analyses of the influence of developmental
mode on phenotypic diversification rates in shorebirds. Proceedings of the Royal Society B 273, 
1619-1624.

Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments
and island effects. Evolution 63, 2017-2030.}
\author{
Gavin Thomas, Rob Freckleton}
\note{
The means are output as treatment contrasts.}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
data(anolis.tree)
data(anolis.data)

## Convert data to class rateData with a rateMatrix object as input
anolis.rateMatrix <- as.rateMatrix(phy=anolis.tree, x="geo_ecomorph", data=anolis.data)

anolis.rateData <- as.rateData(y="Female_SVL", x="geo_ecomorph", 
rateMatrix = anolis.rateMatrix, phy=NULL, data=anolis.data, log.y=TRUE)

## Calculate phylogenetic mean, variance, log likelihood for a model where the first rate category 
#  (small islands) takes position one in the rate vector, group 1 (large island trunk crown and trunk 
#  ground) takes position 2 and so on. The rates in the first and last categories are equal (both 1) 
#  whereas the rate in the second category is twice as fast (2) and the rate in the third category is 
#  ten times slower. Means are allowed to differ. 

# mean only
phyloMean(rateData=anolis.rateData, rate = c(1,2,0.1,1), common.mean = FALSE)

# variance only
phyloVar(rateData=anolis.rateData, rate = c(1,2,0.1,1), common.mean = FALSE) 

# mean, variance and log-likelihood
likRatePhylo(rateData=anolis.rateData, rate = c(1,2,0.1,1), common.mean = FALSE) 


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{optim.likRatePhylo}
\alias{optim.likRatePhylo}
\title{
Maximum likelihood rate estimation for traits and phylogenies 
}
\description{
Function for the maximum likelihood estimation of rate parameters on a trait and phylogeny.}
\usage{
optim.likRatePhylo(rateData, rate = NULL, fixed = NULL, rateMIN = 0.001, rateMAX = 50, common.mean = FALSE, lambda.est=TRUE, meserr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rateData}{
an object of class "rateData"}
  \item{rate}{
a vector of relative rate parameters. The length of the vector is equal to the number of rates being estimated. If \code{rate=NULL} then rates are equal. }
  \item{fixed}{
A vector stating whether each parameter should be allowed to vary (either \code{FALSE} which results in a start value of 1, or a numeric start value) or should be fixed (\code{TRUE}).}
  \item{rateMIN}{
Minimum value for the rate parameters}
  \item{rateMAX}{
Maximum value for the rate parameters}
  \item{common.mean}{
a logical specififying whether each rate category should have its own mean (\code{common.mean=FALSE}) or all categories should have the same mean (\code{common.mean=FALSE}). See Thomas et al. (2009) for a discussion on the impact of assumptions about mean on rate estimates. }
  \item{lambda.est}{
Logical. Fit Pagel's lambda.}

  \item{meserr}{
Logical. Include measurement error. }}

\value{
\item{MLRate }{Maximum likelihood estimates of the rate parameters}
\item{Max.lik }{Maximum (log) likeihood}
\item{AIC }{AIC for maximum likelihood model}
\item{AICc }{AICc for maximum likelihood model}
\item{convergence }{convergence value from \code{optim}}
\item{n.parameters }{Number of parameters in the model (how many means and rate categories)}
}
\references{
Thomas GH, Freckleton RP, & Szekely T. 2006. Comparative analyses of the influence of developmental mode on phenotypic diversification rates in shorebirds. Proceedings of the Royal Society B 273, 1619-1624.

Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments and island effects. Evolution 63, 2017-2030.}
\author{
Gavin Thomas}


%% ~Make other sections like Warning with \section{Warning }{....} ~



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line



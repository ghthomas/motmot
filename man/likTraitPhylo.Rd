\name{likTraitPhylo}
\alias{likTraitPhylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log-likelihood rate estimation for traits and phylogenies
}
\description{
This function calculates the log-likelihood and Brownian (co)variance for a trait(s) 
and a phylogeny using phylogenetically independent contrasts.}
\usage{
likTraitPhylo(y, phy, meserr=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
A matrix of trait data. Rownames must be included and match the taxon names in the phylogeny. Can accept single traits (calculates variance) or multiple traits (calculates variance-covariance matrix). }
  \item{phy}{
An object of class "phylo" (see ape package).}

  \item{meserr}{
A vector (or matrix) of measurement error for each tip. This is only applicable to univariate analyses.}
}
\details{
The "phylo" object must be rooted and fully dichotomous.}
\value{
  \item{brownianVariance }{Brownian variance (or covariance for multiple traits) given the data and phylogeny}
  \item{logLikelihood }{The log-likelihood of the model and data}

}
\references{
Felsenstein J. 1973. Maximum-likelihood estimation of evolutionary trees from continuous characters. Am. J. Hum. Genet. 25, 471-492.

Felsenstein J. 1985. Phylogenies and the comparative method. American Naturalist 125, 1-15.

Freckleton RP & Jetz W. 2009. Space versus phylogeny: disentangling phylogenetic and spatial signals in comparative data. Proc. Roy. Soc. B 276, 21-30. 

}
\author{
Gavin Thomas, Rob Freckleton}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

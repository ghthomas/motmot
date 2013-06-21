\name{phyloCovar}
\alias{phyloCovar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of Brownian (co)variance using independent contrasts.
}
\description{
Calculates the Brownian variance (single trait) or variance-covariance matrix (mutliple traits) using phylogenetically independent contrasts.}
\usage{
phyloCovar(x, phy, estimator = "unbiased")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A continuous trait
}
  \item{phy}{
An object of class "phylo" (see ape package).
}
  \item{estimator}{
Should Brownian variance (or covariance) be based on the unbiased ("unbiased" - default) or maximum likelihood ("ML") estimator.
}
}

\value{

  \item{brownianVariance }{Brownian variance (or covariance for multiple traits) given the data and phylogeny}

}
\references{
Felsenstein J. 1973. Maximum-likelihood estimation of evolutionary trees from continuous characters. Am. J. Hum. Genet. 25, 471-492.

Felsenstein J. 1985. Phylogenies and the comparative method. American Naturalist 125, 1-15.

Freckleton RP & Jetz W. 2009. Space versus phylogeny: disentangling phylogenetic and spatial signals in comparative data. Proc. Roy. Soc. B 276, 21-30. }
\author{
Gavin Thomas, Rob Freckleton
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

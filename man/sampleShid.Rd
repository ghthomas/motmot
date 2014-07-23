\name{sampleShid}
\alias{sampleShid}
\title{
Sample hidden speciation events along branches of a tree
}
\description{
Uses estimated speciation and extinction rates to sample the number of speciation events 'hidden' by subsequent extinction on each branch of a tree following Bokma (2008). For use with the \code{psi} and \code{multipsi} models. 
}
\usage{
sampleShid(phy, la = NULL, mu = NULL, useMean = FALSE)
}
\arguments{
  \item{phy}{
Phylogenetic tree in \code{phylo} format
}
  \item{la}{
Estimate of the rate of speciation "lambda"
}
  \item{mu}{
Estimate of the rate of extinction "mu"
}
  \item{useMean}{
A logical indicating whether to output the average or expected number of hidden speciation events per branch, which may be non-integer (if \code{TRUE}), or to sample an integer number on each branch from a Poisson distribution (if \code{FALSE}, the default)
}
}

\details{
The expected number of hidden speciation events are calculated for each branch given its start and end times, and estimates of lambda and mu which are assumed to be constant across the tree. To properly account for uncertainty in the effect of extinction on the number of nodes affecting each branch of a tree, it may be appropriate to repeat model-fitting on many realizations of \code{Sobs} on the tree of interest (similar to evaluating phylogenetic uncertainty)
}
\value{
Phylogenetic tree in \code{phylo} format, with an added element \code{Sobs}, a vector of numbers of hidden speciation events per branch, in the same order as the branches in the \code{phylo} object
\references{
Bokma, F. 2008. Detection of "punctuated equilibrium" by Bayesian estimation of speciation and extinction rates, ancestral character states, and rates of anagenetic and cladogenetic evolution on a molecular phylogeny. Evolution 62: 2718-2726.
Ingram, T. 2011. Speciation along a depth gradient in a marine adaptive radiation. Proc. R. Soc. B 278: 613-618.
}
\author{
Travis Ingram
}


\keyword{utilities}
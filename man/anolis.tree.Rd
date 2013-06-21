\name{anolis.tree}
\alias{anolis.tree}
\docType{data}
\title{
Anolis phylogeny}
\description{
Phylogeny of Anolis lizards}
\usage{data(anolis.tree)}
\format{
The data are stored as an object of class "phylo". See help pages for the functions \code{read.tree}
and \code{read.nexus} for details of the data structure.
}
\details{
The phylogeny is originally from Nicholson et al. (2005) downloaded from 
http://biosgi.wustl.edu/~lososlab/anolis_mbg_2005/ (this link appears to be broken). The version here was used in Thomas et al. (2009).



}
\source{
Nicholson KE, Glor RE, Kolbe JJ, Larson A, Blair Hedges S, & Losos J. 2005. Mainland colonization
by island lizards. J. Biogeography 32, 929-938.

Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments
and island effects. Evolution 63, 2017-2030.}

\examples{
data(anolis.tree)
plot(anolis.tree)
}
\keyword{datasets}

\name{fairProportions}
\alias{fairProportions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate fair proportions phylogenetic diversity metric}
\description{
Calculate fair proportions phylogenetic diversity metric}
\usage{
fairProportions(phy, nodeCount=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
An object of class "phylo" (see ape package).
}
\item{nodeCount}{Logical - should root to tip node counts be returned (default is FALSE)}
}

\value{
Returns a matrix of fair proportion for all tips in phylogeny and node counts if selected.
}
\references{
Redding, D.W. and Mooers, A.O. (2006). Incorporating evolutionary measures into conservation prioritisation. Conservation Biology, 20, 1670-1678.

Isaac, N.J.B., Turvey, S.T., Collen, B., Waterman, C. and Baillie, J.E.M. (2007). Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE, 2, e296.}
\author{
Gavin Thomas
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

data(anolis.tree)

fp <- fairProportions(anolis.tree)
fpNodes <- fairProportions(anolis.tree, nodeCount=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

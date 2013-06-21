\name{as.rateMatrix}
\alias{as.rateMatrix}
\title{
Conversion among data and phylogeny objects
}
\description{
Function to generate a "rateMatrix" object containing a set of variance covariance matrices. 

Note that \code{as.rateMatrix} calls the CAIC function \code{vcv.array} multiple times and this 
can be slow for large phylogenies (though faster than using the "ape" equivalent \code{vcv.phylo}). }
\usage{
as.rateMatrix(phy, x, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
An object of class "phylo" (see ape package).}
  \item{x}{
The explanatory (discrete) variable used to define the hypothesised rate categories. Can be specified 
as a column number or column name.}
  \item{data}{
A data frame containing (minimally) the x variable as a column with species names as rownames.}
}

\value{

\item{rateMatrix }{An object of class "rateMatrix" - a list of matrices describing the expected 
variances and covariances of between species. Each matrix refers to the variances and covariances 
for a given state of x (see Thomas et al. 2006).}

}
\references{
Thomas GH, Freckleton RP, & Szekely T. 2006. Comparative analyses of the influence of developmental
mode on phenotypic diversification rates in shorebirds. Proceedings of the Royal Society B 273, 
1619-1624.
}
\author{
Gavin Thomas}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
## Read in phylogeny and data from Thomas et al. (2009)
data(anolis.tree)
data(anolis.data)

## Convert data to class rateMatrix
# Not run
# anolis.rateMatrix <- as.rateMatrix(phy=anolis.tree, x="geo_ecomorph", data=anolis.data)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

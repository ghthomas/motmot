\name{make.anc}
\alias{make.anc}
\title{
Create design matrix (internal function)
}
\description{
This is an internal function to generate the design matrix required to define different means for 
each hypothesised rate category.}
\usage{
make.anc(y, x, data = NULL, common.mean = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
The response variable - typically a continuous trait.}
  \item{x}{
The explanatory (discrete) variable used to define the hypothesised rate categories.}
  \item{data}{
A data frame containing (minimally) the x and y variables as columns with species names as rownames.}
  \item{common.mean}{
a logical specififying whether each rate category should have its own mean 
(\code{common.mean=FALSE}) or all categories should have the same mean (\code{common.mean=FALSE}). 
See Thomas et al. (2009) for a discussion on the impact of assumptions about mean on rate estimates.
}
}

\value{

\item{x }{A design matrix}

}
\references{
Thomas GH, Freckleton RP, & Szekely T. 2006. Comparative analyses of the influence of developmental
mode on phenotypic diversification rates in shorebirds. Proceedings of the Royal Society B 273, 
1619-1624.

Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments
and island effects. Evolution 63, 2017-2030.}
\author{
Gavin Thomas}


%% ~Make other sections like Warning with \section{Warning }{....} ~


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{make.likRatePhylo}
\alias{make.likRatePhylo}
\title{
Internal function}
\description{
Internal function. Constructor to allow fixing of rate parameters. }
\usage{
make.likRatePhylo(rateData, fixed, common.mean = FALSE, lambda.est, meserr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rateData}{
an object of class "rateData"}
  \item{fixed}{
A vector stating whether each parameter should be allowed to vary (either \code{FALSE} which 
results in a start value of 1, or a numeric start value) or should be fixed (\code{TRUE}).  }
  \item{common.mean}{
a logical specififying whether each rate category should have its own mean (\code{common.mean=FALSE}) 
or all categories should have the same mean (\code{common.mean=FALSE}). See Thomas et al. (2009) 
for a discussion on the impact of assumptions about mean on rate estimates.  }
  \item{lambda.est}{
Logical - Logical. Fit Pagel's lambda.}
  \item{meserr}{
an object of class "rateData"}
}
\value{
Returns a function to be passed to \code{optim.likRatePhylo}
}

\author{
Gavin Thomas}


%% ~Make other sections like Warning with \section{Warning }{....} ~



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

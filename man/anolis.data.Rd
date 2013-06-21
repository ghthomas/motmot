\name{anolis.data}
\alias{anolis.data}
\docType{data}
\title{
Anolis data}
\description{
Anolis data set from Thomas et al. (2009)}
\usage{data(anolis.data)}
\format{
  A data frame with 164 observations on the following 6 variables.
  \describe{
    \item{\code{Species}}{a factor }
    \item{\code{Island_type}}{a factor with levels \code{LargeIsland} \code{Mainland} 
    \code{SmallIsland}}
    \item{\code{ecomorph}}{a factor with levels \code{CrownGiant} \code{GrassBush} \code{Trunk} 
    \code{TrunkCrown} \code{TrunkGround} \code{Twig} \code{Unique}}
    \item{\code{geo_ecomorph}}{a numeric vector}
    \item{\code{Female_SVL}}{a numeric vector of female snout vent length}
    \item{\code{Male_SVL}}{a numeric vector of male snout vent length}
  }
}

\source{
Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments 
and island effects. Evolution 63, 2017-2030.}

\examples{
data(anolis.data)
}
\keyword{datasets}
\name{pic.motmot}
\alias{pic.motmot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Phylogenetically independent contrasts (internal)}
\description{
Calculates phylogenetically independent contrasts.  }
\usage{
pic.motmot(x, phy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A matrix of trait values with taxon names as rownames.
}
  \item{phy}{
An object of class "phylo" (see ape package).
}

}
\details{
Extracts values for contrasts, expected variances using contrasts by calling \code{pic}. }
\value{
  \item{contr }{A matrix with two columns containing raw contrasts in the first column and their expected variances in the second column.}
  \item{root.v }{Expected variances of branches either side of root}
  \item{V }{Expected variance at the root }

}
\references{
Felsenstein J. 1985. Phylogenies and the comparative method. American Naturalist 125, 1-15.
}
\author{
Gavin Thomas, Rob Freckleton, Emmanuel Paradis}


%% ~Make other sections like Warning with \section{Warning }{....} ~



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

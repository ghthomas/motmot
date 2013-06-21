\name{node.descendents}
\alias{node.descendents}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify nodes and tips descended from a node (internal function).
}
\description{
Obtains a vector of the tips and nodes subtending from a node in a phylogeny.}
\usage{
node.descendents(x, phy, tip.labels=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A positive integer 
}
  \item{phy}{
An object of class "phylo" (see ape package).
}

  \item{tip.labels}{
Logical - output tip.labels?
}

}
\details{
This function is stolen from the clade.members function in the CAIC package but returns both node and tip id's.}
\value{
Returns a vector of node and tip ids descended from the tip(s) "x". If tip.labels=TRUE then returns a list of node ids and tip labels.
}

\author{
Gavin Thomas, David Orme}


%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

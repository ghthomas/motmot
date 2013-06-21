\name{cladeIdentity}
\alias{cladeIdentity}
\title{
Identify branches (including tips) descended from a node (internal function).
}
\description{
Internal function to get presence absence of descendent branches from a vector of node numbers. The descendents include the branch leading to the focal node (i.e. node defines the stem group no crown group)}
\usage{
cladeIdentity(phy, nodeIDs, cladeMembersObj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
An object of class "phylo" (see ape package).
}
  \item{nodeIDs}{
Vector of node numbers (positive integers).}
  \item{cladeMembersObj}{
Matrix of clade membership}
}
\details{
The function returns a matrix of unique presences given the selected node. If the selected nodes are nested then presences are only recorded for the least inclusive node.
}
\value{

\item{matrix }{Matrix of unique presences for each node id}

}

\author{
Gavin Thomas}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

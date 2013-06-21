\name{traitMedusaSummary}
\alias{traitMedusaSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify shifts in the rate of trait diversification
}
\description{
Summarises phenotypic rate variation on phylogenies. }
\usage{
traitMedusaSummary(traitMedusaObject = NULL, cutoff = 4, AICc = TRUE, lowerBound=1e-7, upperBound=200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{traitMedusaObject}{
Output of a medusa analysis in transformPhylo.ML }
  \item{cutoff}{
Cutoff value for differences in AIC scores when comparing models. More complex models with an AIC score more than this number of units lower than simpler models are retained (as per runMedusa in the geiger package).}
  \item{AICc}{
If true, AICc is used instead of AIC.}
  \item{lowerBound}{
Minimum value for parameter estimates.}
  \item{upperBound}{
Maximum value for parameter estimates. }

}
\details{
This functions summarises the output of a "medusa" model in transformPhylo.ML (see below). The best overall model is chosen based on AIC (or AICc if AICc=TRUE). The cut-off point for improvement in AIC score between successively more complex models can be defined using cutoff. The default cutoff is 4 but this is somewhat arbitrary and a "good" cut-off may well vary between data sets so it may well be worth exploring different cutoffs.

Summarises fits of "medusa" models ("clade" models generated without any a priori assertion of the location of phenotypic diversification rate shifts). It uses the same approach as the runMedusa function in the geiger package (runMedusa tests for shifts in the rate of lineage diversification). The algorithm first fits a constant-rate Brownian model to the data, it then works iteratively through the phylogeny fitting a two-rate model at each node in turn. Each two-rate model is compared to the constant rate model and the best two-rate model is retained. Keeping the location of this rate shift intact, it then repeats the procedure for a three-rate model and so on. The maximum number of rate shifts can be specified a priori using nSplits. Limits can also be applied to the size (species richness) of clades on which to infer new rate shifts using minCladeSize. This can be useful to enable large trees to be handled but should be used cautiously since specifiying a large minimum clade size may result in biologically interesting nested rate shifts being missed. Equally, very small clade sizes may provide poor estimates of rate that may not be informative.}
\value{
    \item{ModelFit }{Summary of the best optimal rate shift model.}
    \item{Rates }{Summary of the rate parameters from the best rate shift model.}
    \item{optimalTree }{A phylo object with branch lengths scaled relative to rate.}


}
\references{
Alfaro ME, Santini F, Brock CD, Alamillo H, Dornburg A, Carnevale G, Rabosky D & Harmon LJ. 2009. Nine exceptional radiations plus high turnover explain species diversity in jawed vertebrates. PNAS 106, 13410-13414.

O'Meara BC, Ane C, Sanderson MJ & Wainwright PC. 2006. Testing for different rates of continuous trait evolution using likelihood. Evolution 60, 922-933

Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments and island effects. Evolution 63, 2017-2030.}
\author{
Gavin Thomas
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
data(anolis.tree)
data(anolis.data)

# anolis.data is not matrix and contains missing data so put together matrix of # relevant traits (here female and male snout vent lengths) and remove species 
# with missing data from the matrix and phylogeny
anolisSVL <- data.matrix(anolis.data)[,c(5,6)]
anolisSVL[,1] <- log(anolisSVL[,1])
anolisSVL[,2] <- log(anolisSVL[,2])

tree <- drop.tip(anolis.tree, names(attr(na.omit(anolisSVL), "na.action")))
anolisSVL <- na.omit(anolisSVL)

# Identify rate shifts and print and plot results
# Not run
#anolisSVL_MEDUSA <- transformPhylo.ML(anolisSVL, phy=tree, model="tm1", minCladeSize=20, nSplits=3)

#anolisSVL_MEDUSA_out <- traitMedusaSummary(anolisSVL_MEDUSA)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

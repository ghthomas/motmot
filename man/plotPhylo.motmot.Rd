\name{plotPhylo.motmot}
\alias{plotPhylo.motmot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tree plotting for rates}
\description{
Plots trees with colours based on rates of trait evolution. Also provides simple coloured plotting for trait values using the "ace" function in the ape library.}
\usage{
plotPhylo.motmot(phy, x=NULL, traitMedusaObject=NULL, reconType="rates", type = "phylogram", use.edge.length = TRUE, show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black", edge.width = 1, edge.lty = 1, font = 3, cex = par("cex"), adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE, label.offset = 0.5, underscore = FALSE, x.lim = NULL, y.lim = NULL, direction = "rightwards", lab4ut = "horizontal", tip.color = "black", palette="hotspot.colors")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
	An object of class "phylo" (see ape package).}
  \item{x}{
	A matrix of trait values.}
  \item{traitMedusaObject}{
	Output from traitMedusaSummary.}
  \item{reconType}{
	Colour branches according to rate shifts ("rates" - requires traitMedusaObject) or ancestral state reconstruction ("picReconstruction"  - requires x). }
  \item{palette}{
	Defines the colour scheme with four options: hotspot.colors (red to blue), heat.colors (yellow to red), cool.colors (blues), combi.colors (yellows to reds 			and blues)}
  \item{type}{a character string specifying the type of phylogeny to be
    drawn; it must be one of "phylogram" (the default), "cladogram",
    "fan", "unrooted", "radial" or any unambiguous abbreviation of
    these.}
  \item{use.edge.length}{a logical indicating whether to use the edge
    lengths of the phylogeny to draw the branches (the default) or not
    (if \code{FALSE}). This option has no effect if the object of class
    \code{"phylo"} has no `edge.length' element.}
  \item{show.tip.label}{a logical indicating whether to show the tip
    labels on the phylogeny (defaults to \code{TRUE}, i.e. the labels
    are shown).}
  \item{show.node.label}{a logical indicating whether to show the node
    labels on the phylogeny (defaults to \code{FALSE}, i.e. the labels
    are not shown).}
  \item{edge.color}{a vector of mode character giving the colours used
    to draw the branches of the plotted phylogeny. These are taken to be
    in the same order than the component \code{edge} of \code{phy}. If
    fewer colours are given than the length of \code{edge}, then the
    colours are recycled.}
  \item{edge.width}{a numeric vector giving the width of the branches of
    the plotted phylogeny. These are taken to be in the same order than
    the component \code{edge} of \code{phy}. If fewer widths are given
    than the length of \code{edge}, then these are recycled.}
  \item{edge.lty}{same than the previous argument but for line types;
    1: plain, 2: dashed, 3: dotted, 4: dotdash, 5: longdash, 6: twodash.}
  \item{font}{an integer specifying the type of font for the labels: 1
    (plain text), 2 (bold), 3 (italic, the default), or 4 (bold
    italic).}
  \item{cex}{a numeric value giving the factor scaling of the tip and
    node labels (Character EXpansion). The default is to take the
    current value from the graphical parameters.}
  \item{adj}{a numeric specifying the justification of the text strings
    of the labels: 0 (left-justification), 0.5 (centering), or 1
    (right-justification). This option has no effect if \code{type =
      "unrooted"}. If \code{NULL} (the default) the value is set with
    respect of \code{direction} (see details).}
  \item{srt}{a numeric giving how much the labels are rotated in degrees
    (negative values are allowed resulting in clock-like rotation); the
    value has an effect respectively to the value of
    \code{direction} (see Examples). This option has no effect if
    \code{type = "unrooted"}.}
  \item{no.margin}{a logical. If \code{TRUE}, the margins are set to
    zero and the plot uses all the space of the device (note that this
    was the behaviour of \code{plot.phylo} up to version 0.2-1 of `ape'
    with no way to modify it by the user, at least easily).}
  \item{root.edge}{a logical indicating whether to draw the root edge
    (defaults to FALSE); this has no effect if `use.edge.length = FALSE'
    or if `type = "unrooted"'.}
  \item{label.offset}{a numeric giving the space between the nodes and
    the tips of the phylogeny and their corresponding labels. This
    option has no effect if \code{type = "unrooted"}.}
  \item{underscore}{a logical specifying whether the underscores in tip
    labels should be written as spaces (the default) or left as are (if
    \code{TRUE}).}
  \item{x.lim}{a numeric vector of length one or two giving the limit(s)
    of the x-axis. If \code{NULL}, this is computed with respect to
    various parameters such as the string lengths of the labels and the
    branch lengths. If a single value is given, this is taken as the
    upper limit.}
  \item{y.lim}{same than above for the y-axis.}
  \item{direction}{a character string specifying the direction of the
    tree. Four values are possible: "rightwards" (the default),
    "leftwards", "upwards", and "downwards".}
  \item{lab4ut}{(= labels for unrooted trees) a character string
    specifying the display of tip labels for unrooted trees: either
    \code{"horizontal"} where all labels are horizontal (the default),
    or \code{"axial"} where the labels are displayed in the axis of the
    corresponding terminal branches. This option has an effect only if
    \code{type = "unrooted"}.}
  \item{tip.color}{the colours used for the tip labels, eventually
    recycled (see examples).} 

}

\value{

Returns a data frame of colours used in plot along with rate (or ancestral state) range for each colour. 

}

\author{
Gavin Thomas}

\examples{
# Data and phylogeny
data(anolis.tree)
data(anolis.data)

# anolis.data is not matrix and contains missing data so put together matrix of # relevant traits (here female and male snout vent lengths) and remove species 
# with missing data from the matrix and phylogeny
anolisSVL <- data.matrix(anolis.data)[,c(5,6)]
anolisSVL[,1] <- log(anolisSVL[,1])
anolisSVL[,2] <- log(anolisSVL[,2])

tree <- drop.tip(anolis.tree, names(attr(na.omit(anolisSVL), "na.action")))
anolisSVL <- na.omit(anolisSVL)

# Identify rate shifts and print and plot results with upto three rate shifts and minimum clade size of 20.
# Not run
# anolisSVL_MEDUSA <- transformPhylo.ML(anolisSVL, phy=tree, model="tm1", minCladeSize=10, nSplits=2)

# anolisSVL_MEDUSA_out <- traitMedusaSummary(anolisSVL_MEDUSA, cutoff=3, AICc=FALSE)

# Plot rate shifts
# colours <- plotPhylo.motmot(phy=tree, traitMedusaObject = anolisSVL_MEDUSA_out,  reconType = "rates", type = "fan", cex=0.6, edge.width=3)

# Plot pic ancestral state reconstruction for female SVL
# colours <- plotPhylo.motmot(phy=tree, x=anolisSVL[,1], palette="hotspot.colors", edge.width=4, cex=0.8, show.tip.label=TRUE, adj=0.5, label.offset=2, reconType="picReconstruction")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

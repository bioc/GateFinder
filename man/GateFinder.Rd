\name{GateFinder}
\alias{GateFinder}
\title{
GateFinder
}
\description{
Given a vector of cluster memberships for a cell population, identifies a sequence of gates (polygon filters on 2D scatter plots) for isolation of that cell type.
}
\usage{
GateFinder(x, targetpop, update.gates=FALSE, max.iter=2, beta=1, outlier.percentile=0.05, subsample=length(targetpop), nstart=1, update.org.data=TRUE, randomize=(nstart>1), selection.criteria='best', unimodalitytest=TRUE, predimx=NULL, predimy=NULL, convex=TRUE, alpha=5)
}

\arguments{
  \item{x}{
     A flowFrame or an expression matrix in which columns are markers and rows are cells.
}
  \item{targetpop}{
     A vector of logical values one for each cell (TRUE = cells in the target population). If instead a vector of integers are supplied each integer value will be treated as a separate celltype of interest and a list of gating strategies will be returned.
}
  \item{update.gates}{
     A boolean value indicating if the polygon gates should be updated after each gating step. update.gates=TRUE makes the analysis slower.
}
  \item{max.iter}{
  The number of requested gating steps.
}
  \item{beta}{
  A positive real value which control the trade-off between precision and recall in the F-measure calculation. Values smaller than 1 (and closer to 0) emphasize recall values and values larger than 1 emphasize precision.
}
  \item{outlier.percentile}{
  The percentile of the empirical distribution of each 2d distribution (each scatter plot) that should be excluded before calculating the polygon gate. If a vector, all provided numbers will be tested and the one with the highest F-measure will be used.
}
  \item{subsample}{
  The number of cells (integer) that should be randomly selected for calculation of the gating strategy (the gating strategy can be applied to all cells therefore the resulting populations will be based on all cells (see \code{update.org.data})).
}
  \item{nstart}{
  The number of randomized runs (integer). The results from the best (or median) randomized run will be used. See \code{selection.criteria}).
}
  \item{update.org.data}{
  A boolean controling weather the gating strategy calculated using a random subset of the cells should be applied to all cells or not.
}
  \item{randomize}{
  A boolean value to control if selection of the gates should be randomized. If TRUE the selection probability of each gate will be proportional to it's fmeasure. Otherwise the best gate will be selected for each step.
}
  \item{selection.criteria}{
  A string with values of either 'best' or 'median'. This determines if the run with the best or median fmeasure should be used as the final gating strategy.
}
  \item{unimodalitytest}{
  A boolean value. If TRUE the unimodality of the first principal component will be tested using a dip test and a warning is issued for p < 0.05.
}
  \item{predimx}{
A vector of marker numbers for the x-axes of a pre-determined gating strategy.
}
  \item{predimy}{
A vector of marker numbers for the y-axes of a pre-determined gating strategy.
}
  \item{convex}{
A boolean value indicating if the target population is expected to be convex (for outlier removal purposes).
}

\item{alpha}{
	alpha-hull threshold for non-convex gates.
}

}
\value{
\item{GatingProjection}{A GatingProjection object.}
}
\references{
Filzmoser, Peter, Ricardo Maronna, and Mark Werner. "Outlier identification in high dimensions." Computational Statistics & Data Analysis 52, no. 3 (2008): 1694-1711.
}
\author{
Nima Aghaeepour <naghaeep@gmail.com> and Erin F. Simonds <erin.simonds@gmail.com>.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{


library(flowCore)

data(LPSData) 
##Select the target population. In this case cells with those with a pP38 expression (dimension 34) of higher than 3.5.
targetpop <- (exprs(rawdata)[,34]>3.5)

##Subset the markers that should be considered for gating.
x=exprs(rawdata)[,prop.markers]
colnames(x)=marker.names[prop.markers]

##Run GateFinder.
ans=GateFinder(x, targetpop)

##Make the plots.
plot (x, ans, c(2,3), targetpop)
plot(ans)


##Alternatively, using a flowFrame:
x=new('flowFrame', exprs=x)
ans=GateFinder(x, targetpop)

##Now you can use the gates and filters to subset the flowFrame. E.g.:
split(x, ans@flowEnv$Filter2)

##This function relies on an EXPERIMENTAL feature in flowUtils. Please be cautious when replying on this.
##Don't run without the optional flowUtils package installed.
##To write the gates into a GatingML file:
##library(flowUtils)
##write.gatingML(ans@flowEnv, 'GatingML.xml')

}

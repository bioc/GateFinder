\name{plot.GateFinder}
\alias{plot.GateFinder}
\title{
plot.GateFinder
}
\description{
Creates a scatter plot for each of the gating steps.
}
\usage{
plot.GateFinder(x, y, ncolrow=c(1,max(targetpop)), targetpop=NULL, beta=NULL,  cexs=NULL, cols=NULL, subsample=length(targetpop), max.iter=length(y@gates), pot=TRUE, xlim=NULL, ylim=NULL, asinh.axis=FALSE,...)
}

\arguments{
  \item{x}{
     A flowFrame or an expression matrix in which columns are markers and rows are cells.
}
  \item{y}{
     A GatingProjection object.
}
  \item{ncolrow}{
     A vector of length 2 indicating the desired number of rows and columns in the plot.
}
  \item{targetpop}{
  The target cell type.
}
  \item{beta}{
  A positive real value which control the trade-off between precision and recall in the F-measure calculation. Values smaller than 1 (and closer to 0) emphasize recall values and values larger than 1 emphasize precision.
}
  \item{cexs}{
  A vector of length 3 indicating the point sizes for 1-previously excluded cells 2-non-selected cells 3-selected cells.
}
  \item{cols}{
  A vector of length 3 indicating the point colors for 1-previously excluded cells 2-non-selected cells 3-selected cells.
}
  \item{subsample}{
  The number of randomized runs (integer). The results from the best (or median) randomized run will be used. See \code{selection.criteria}).
}
  \item{max.iter}{
  A boolean controling weather the gating strategy calculated using a random subset of the cells should be applied to all cells or not.
}
  \item{pot}{
  A boolean value. If true, the points of interest will be ploted on top of other points to increase visibility.
}
  \item{xlim}{
  Static x-axis limits for the plots (vector of length 2).
}
  \item{ylim}{
  Static y-axis limits for the plots (vector of length 2).
}
  \item{asinh.axis}{
A boolean value indicating if asinh axis ticks should be plotted (usually used for mass cytometry data).
}
  \item{...}{
Other arguments passed to the plot function.
}

}

\value{
\item{Plot}{A GateFinder plot.}
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

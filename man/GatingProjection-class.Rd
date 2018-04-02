\name{GatingProjection-class}
\Rdversion{1.1}
\docType{class}
\alias{GatingProjection}
\alias{GatingProjection-class}
\alias{plot,GatingProjection,ANY-method}
\alias{plot,matrix,GatingProjection-method}

\title{Class \code{"GatingProjection"}}
\description{
An object that stores the final gating projections as well as the scores calculated for each step.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GatingProjection", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{fmeasure}:}{A vector of F-measure values for each step of the identified gating strategy.}
    \item{\code{precision}:}{A vector of precision values for each step of the identified gating strategy.}
    \item{\code{recall}:}{A vector of recall values for each step of the identified gating strategy.}
    \item{\code{dimx}:}{A vector of marker indexes for the x-axis of each step of the identified gating hierarchy.}
    \item{\code{dimy}:}{A vector of marker indexes for the y-axis of each step of the identified gating hierarchy.}
    \item{\code{gates}:}{A list of polygon gates for each step of the identified gating hierarchy.}
    \item{\code{pops}:}{A list of vectors representing the cell population memberships for each step of the identified hierarchy.}
    \item{\code{subsampleindex}:}{A vector of the indexes of the selected subsample of cells (if applicable).}
    \item{\code{fmeasures}:}{A vector of F-measure values of multiple randomized attempts (if applicable).}
    \item{\code{flowEnv}:}{An envirnoment for flowCore's polygon gates and intersect filters.}
  }
}

\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "GatingProjection", y = "ANY")}:
Plot of F-measure, precision, and recall values of each gating step. }
    \item{plot}{\code{signature(x = "matrix", y = "GatingProjection")}:
Scatter plots of the raw data (from matrix x) for each step of the gating strategy. Gray dots represent cells that were removed in the previous step. Red dots represent the target cells.
}
	 }
}

\author{
Nima Aghaeepour <naghaeep@gmail.com> and Erin F. Simonds <erin.simonds@gmail.com>.
}



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

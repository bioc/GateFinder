\name{GateFinder-package}
\alias{GateFinder-package}
%\alias{GateFinder}
\docType{package}
\title{
GateFinder
}
\description{
Given a vector of cluster memberships for a cell population, identifies a sequence of gates (polygon filters on 2D scatter plots) for isolation of that cell type.
}
\details{
\tabular{ll}{
Package: \tab GateFinder\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-12-21\cr
License: \tab Artistic-2.0\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Nima Aghaeepour <naghaeep@gmail.com> and Erin F. Simonds <erin.simonds@gmail.com>
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

\name{LPSData}
\alias{LPSData}
\alias{rawdata}
\alias{marker.names}
\alias{prop.markers}
\docType{data}
\title{LPSData}
\description{
  A dataset of two sets of scores (particularly, correlation with protection against HIV and overlap with the Naive T-cell population) assigned to immunophenotypes measured by flow cytometry. 10 markers were measured: KI-67, CD28, CD45RO, CD8, CD4, CD57, CCR5, CD27, CCR7, and CD127.
}
\usage{data(LPSData)}
\details{
  This dataset consists of a matrix and two vectors:

 \code{rawdata} The transformed expression values extracted from the original FCS file.

  \code{prop.markers} the indexes of markers that should be considered for the gating strategy.

\code{marker.names} name of all markers (columns of matrix rawdata).


}
%\source{
%}

\author{
  Nima Aghaeepour \email{<naghaeep@gmail.com>} and Erin F. Simonds \email{<erin.simonds@gmail.com>}.
}

\references{
Bendall, Sean C., et al. "Single-cell mass cytometry of differential immune and drug responses across a human hematopoietic continuum." Science 332.6030 (2011): 687-696.
}

\examples{

library(flowCore)

data(LPSData)

plot(exprs(rawdata)[,15:16],xlab=marker.names[15],ylab=marker.names[16])
}
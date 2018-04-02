test.GateFinder <- function(){
    library(RUnit)
    library(GateFinder)
    library(flowCore)
    dimx=c(6, 3, 1, 7, 5)[1:2]
    dimy=c(8, 12, 4, 11, 9)[1:2]
    fmeasure=c(0.8948949, 0.9006211, 0.8895899, 0.8789809, 0.8636364)[1:2]
    precision=c(0.9197531, 0.9602649, 0.9657534, 0.9650350, 0.9708029)[1:2]
    recall=c(0.8713450, 0.8479532, 0.8245614, 0.8070175, 0.7777778)[1:2]

    data(LPSData)
    targetpop <- (exprs(rawdata)[,34]>3.5)
    x=exprs(rawdata)[,prop.markers]
    colnames(x)=marker.names[prop.markers]
    ans=GateFinder(x, targetpop)

    checkEquals(as.vector(ans@dimx),dimx)
    checkEquals(as.vector(ans@dimy),dimy)
    checkEquals(ans@fmeasure,fmeasure,tolerance=10^-6)
    checkEquals(ans@precision,precision,tolerance=10^-6)
    checkEquals(ans@recall,recall,tolerance=10^-6)
}



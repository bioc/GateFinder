##Main function
GateFinder <- function(x, targetpop, update.gates=FALSE, max.iter=2, beta=1, outlier.percentile=0.05, subsample=length(targetpop), nstart=1, update.org.data=TRUE, randomize=(nstart>1), selection.criteria='best', unimodalitytest=TRUE, predimx=NULL, predimy=NULL, convex=TRUE, alpha=5){
  
  ##multiple target celltypes are handled separately (and recursively) here.
  if (max(targetpop)>1){
    ##if (FALSE){
    ##bootstrapping subspaces
    ##browser()
    prop.markers=seq(ncol(x))
    mat=matrix(0,length(prop.markers),100)
    for (i in 1:ncol(mat)){
      subsampleindex <- sample(seq(nrow(x)),subsample)
      pm <- sample(seq(length(prop.markers)),4)
      temp=.scoreSubSpace(x[subsampleindex,prop.markers[pm]],targetpop[subsampleindex], beta=beta)
      mat[pm,i]=temp
    }
    ##rownames(mat)=marker.names[prop.markers]
    ##variable importance extracted from the bootstraps
    imp=rowMeans(mat)
    ##prop.markers2 are the markers selected for the gating strategies.
    prop.markers2=prop.markers[order(abs(imp),decreasing=TRUE)][1:round(4*2)]
    orgimpmarkers=prop.markers2
    scores=list()
    markers=list()
    while (length(prop.markers2)>max.iter*2){
      tempscores=vector()
      for (i in 1:length(prop.markers2)){
        temp=prop.markers2[-i]
        tempscores[i]=.scoreSubSpace(x[subsampleindex,temp],targetpop[subsampleindex], beta=beta)
      }
      prop.markers2 <- prop.markers2[-which.max(tempscores)]
      markers[[length(markers)+1]]=prop.markers2
      scores[[length(scores)+1]]=tempscores
    }
    ##subsample the data
    x=x[,prop.markers2]
    ##create and return the gating strategies
    gas=list()
    for (c in seq(max(targetpop))){
      localtargetpop <- (targetpop==c)
      gas[[c]]=GateFinder(x, localtargetpop, update.gates=update.gates, max.iter=max.iter, beta=beta, outlier.percentile=outlier.percentile, subsample=subsample, nstart=nstart, update.org.data=update.org.data, randomize=randomize, selection.criteria=selection.criteria, unimodalitytest=unimodalitytest, predimx=predimx, predimy=predimy, convex=convex, alpha=alpha)
    }
    return (gas)
  }
  
  if (!is.null(predimx))
    update.gates=TRUE
  
  if (class(x)=='flowFrame')
    x=exprs(x)
  
  if (unimodalitytest){
    pc <- prcomp(x[targetpop,])
    pv=dip.test(pc$x[,1])$p.value
    if(pv<0.05)
      warning(sprintf('The target cell population does not appear to be single modal. P-value: %f',pv))
  }
  
  ##bootsrapping
  if (nstart>1){
    ga=list()
    fmeasures=vector()
    for (i in 1:nstart){
      if (subsample==length(targetpop))
        stop('For nstart > 1 you need to set the subsample parameter to a value less that the total number of cells.')
      ga[[i]]=GateFinder(x, targetpop, update.gates=update.gates, max.iter=max.iter, beta=beta, outlier.percentile=outlier.percentile, subsample=subsample, nstart=1, update.org.data=TRUE,predimx=predimx, predimy=predimy)
      fmeasures[i]=mean(ga[[i]]@fmeasure)
    }
    ga[[which.max(fmeasures)]]@fmeasures=fmeasures
    if (selection.criteria=='best')
      return(ga[[which.max(fmeasures)]])
    if (selection.criteria=='median')
      return(ga[[which(fmeasures==median(fmeasures))]])            
  }
  
  ##subsampling
  org.x=x
  org.targetpop=targetpop
  subsampleindex=NA
  if (subsample < length(targetpop)){
    subsampleindex <- sample(seq(nrow(x)),subsample)
    x=x[subsampleindex,]
    targetpop=targetpop[subsampleindex]
  }
  if (subsample > length(targetpop)){
    subsampleindex <- sample(seq(nrow(x)),subsample-length(targetpop))
    x=rbind(x,x[subsampleindex,])
    targetpop=c(targetpop,targetpop[subsampleindex])
  }
  
  ##gating
  gates=NULL
  fmeasure=vector()
  recall=vector()
  precision=vector()
  dimx=vector()
  dimy=vector()
  currentpop=rep(TRUE,length(targetpop))    
  finalgates=list()
  pops=list()
  for (i in 1:max.iter){
    if (i==1 | update.gates)
      gates=.getGates(x[currentpop,],targetpop[currentpop], skipdims=c(dimx,dimy),outlier.percentile=outlier.percentile, beta, predimx=predimx[i], predimy=predimy[i], convex=convex, alpha=alpha)
    temp=.getScores(x, targetpop, gates=gates, skipdims=c(dimx,dimy), currentpop=currentpop, beta, randomize=randomize, predimx=predimx[i], predimy=predimy[i])
    if (is.na(temp)){
      if (i==1)
        return(NULL)
      break
    }
    fmeasure[i]=temp$fmeasure
    recall[i]=temp$recall
    precision[i]=temp$precision
    dimx[i]=temp$i
    dimy[i]=temp$j
    finalgates[[i]]=gates[[temp$i]][[temp$j]]
    currentpop=temp$newpop
    pops[[i]]=currentpop
  }
  
  flowEnv <- new.env()
  currentfilter=list()
  for (i in 1:length(finalgates)){
    sqrcut <- finalgates[[i]]
    colnames(sqrcut) <- colnames(x)[c(dimx[i],dimy[i])]
    flowEnv[[sprintf('Gate%d', i)]] <- polygonGate(filterId=sprintf('Gate%d', i), .gate= sqrcut)
    currentfilter[[i]] <- flowEnv[[sprintf('Gate%d', i)]]
    andGate <- new("intersectFilter", filterId=sprintf('Filter%d', i), filters=currentfilter)
    flowEnv[[sprintf('Filter%d', i)]] <- andGate
  }
  
  names(dimx)=colnames(x)[dimx]
  names(dimy)=colnames(x)[dimy]
  
  ga=new('GatingProjection',fmeasure=fmeasure, precision=precision, recall=recall, dimx=dimx, dimy=dimy, gates=finalgates, pops=pops, subsampleindex=subsampleindex, flowEnv=flowEnv)
  
  ##upsampling
  if((subsample<length(org.targetpop)) & update.org.data){
    ga=.gapply(org.x, ga, org.targetpop, beta=beta)
  }
  return(ga)
}

##Given a cell population (targetpop) and a matrix (x), construct polygon gates for all pairs of dimensions (except those in skipdims) using a convex hull of all cells except 0.05 outliers.
.getGates <- function(x, targetpop, skipdims=NA, outlier.percentile, beta ,predimx=NULL, predimy=NULL, convex=TRUE, alpha=5){
  n=dim(x)[2]
  gates=list()
  for (i in 1:(n-1)){
    if (i %in% skipdims)
      next
    gates[[i]]=list()
    for (j in (i+1):n){
      gates[[i]][[j]]=list()
      if (!is.null(predimx))
        if (i!=predimx)
          next
      if (!is.null(predimy))
        if (j!=predimy)
          next
      pts=x[targetpop,c(i,j)]
      pi=NULL
      ans=rep(1,nrow(pts))
      tryCatch({
        if(convex)
          ans=pcout(pts)$wscat
        if(!convex)
          ans=.flowFPOutliers(x[,c(i,j)], targetpop)[targetpop]
      },  error = function(e){
        pi=1:(dim(pts)[1])
        gates[[i]][[j]]=pts[pi[chull(pts[pi,])],]                    
        warning(sprintf('Artificial Distribution Detected in Scatter Plot %s vs %s. Outlier detection was disabled for this plot. If you end up getting an error message please consider using the function jitter() or something similar.',colnames(x)[i], colnames(x)[j]))
      })
      if (is.null(pi)){
        ##pi=which(ans>outlier.percentile)
        ##gates[[i]][[j]]=pts[pi[chull(pts[pi,])],]
        gates[[i]][[j]]=.getBestPoly(x, i, j, ans, targetpop, beta, outlier.percentile, convex=convex, alpha=alpha)
      }
    }
  }
  return(gates)
}

##outlier scoring using flowFP's probability binning
.flowFPOutliers <- function(x, targetpop){
  fp<-flowFPModel(new('flowFrame', exprs=x[targetpop,]))
  fpa<-flowFP(new('flowFrame', exprs=x[,]), fp)
  fpt<-flowFP(new('flowFrame', exprs=x[targetpop,]), fp)
  fpnt<-flowFP(new('flowFrame', exprs=x[!targetpop,]), fp)
  ##bins=which(fpt@counts/fpnt@counts>0.1)
  scores=(fpt@counts/fpnt@counts)
  ans=(scores[fpa@tags[[1]]])
  return(ans)        
}

##spline fitting for smooth gate boundaries
.spline.poly <- function(xy, vertices, k=3, ...) {
  # Assert: xy is an n by 2 matrix with n >= k.
  # Wrap k vertices around each end.
  n <- dim(xy)[1]
  if (k >= 1) {
    data <- rbind(xy[(n-k+1):n,], xy, xy[1:k, ])
  } else {
    data <- xy
  }
  # Spline the x and y coordinates.
  data.spline <- spline(1:(n+2*k), data[,1], n=vertices, ...)
  x <- data.spline$x
  x1 <- data.spline$y
  x2 <- spline(1:(n+2*k), data[,2], n=vertices, ...)$y
  # Retain only the middle part.
  cbind(x1, x2)[k < x & x <= n+k, ]
}


.binsearch <- function(){
  b=0
  e=1
  sb=test(b)
  se=test(e)
  c=(b+e)/2
  ce=test(c)    
  while(min(abs(c-b), abs(c-e))<0.05){
    
  }
}

##given a range of percentiles, find the polygon that best fits the data in the two respective dimensions.
.getBestPoly <- function(x, dimx, dimy, scores ,targetpop, beta, percentiles, convex=TRUE, alpha=5){
  orgalpha=alpha
  pts=x[targetpop,c(dimx,dimy)]
  precision=vector()
  recall=vector()
  fmeasure=vector()
  if (is.null(scores)){
    stop('Outlier scores are missing')
  }
  gate=list()
  for (i in 1:length(percentiles)){
    pi=which(scores>quantile(scores, percentiles[i]))
    pi=pi[!duplicated(pts[pi,])]
    gate[[i]]=NA
    if(length(pi)<3){
      fmeasure=0
      next
    }
    if (convex)
      gate[[i]]=pts[pi[chull(pts[pi,])],]
    if (!convex){
      b=NULL
      ahull.obj=NULL
      oldalpha=alpha
      tryCatch({
        tpts=pts[pi,]
        tpts=tpts[sample(seq(nrow(tpts)),min(nrow(tpts), 500)),]
        ahull.obj=ahull(tpts, alpha=alpha)
        b=.ah2sp(ahull.obj)
        num=length(b[1]@polygons[[1]]@Polygons)
        cnt=1
        while(num>1){
          alpha=alpha*2
          ahull.obj=ahull(tpts, alpha=alpha)
          b=.ah2sp(ahull.obj)
          oldnum=num
          num=length(b[1]@polygons[[1]]@Polygons)
          if (num>oldnum | alpha>500){
            gate[[i]]=pts[pi[chull(pts[pi,])],]
            alpha=orgalpha
            break
          }
          cnt=cnt+1
        }
      },  error = function(e){
        warning(sprintf('Alpha-hull construction in scatter plot %s vs %s could not be constructed. Adding a small amount of noise to try again.',colnames(x)[dimx], colnames(x)[dimy]))
        gate[[i]]=NA
        ahull.obj=ahull(jitter(tpts), alpha=alpha)
        b<<-.ah2sp(ahull.obj)
        num=length(b[1]@polygons[[1]]@Polygons)
        cnt=1
        while(num>1){
          alpha=alpha*2
          ahull.obj=ahull(jitter(tpts), alpha=alpha)
          b<<-.ah2sp(ahull.obj)
          oldnum=num
          num=length(b[1]@polygons[[1]]@Polygons)
          if (num>oldnum | alpha>500){
            gate[[i]]=pts[pi[chull(pts[pi,])],]
            alpha=orgalpha
            break
          }
          cnt=cnt+1
        }
      })
      
      if (!is.na(gate[[i]]))
        next
      ##bsx=ahull.obj$xahull[(as.vector(ahull.obj$arcs[,7])), ]
      bsx=b[1]@polygons[[1]]@Polygons[[1]]@coords
      if(nrow(bsx)>=3)
        gate[[i]]=.spline.poly(bsx,100)
    }
    foundpop=as.logical(.inpoly(x[,dimx],x[,dimy],list(x=gate[[i]][,1],y=gate[[i]][,2])))
    tp=length(which(targetpop & foundpop))
    fp=length(which(!targetpop & foundpop))
    fn=length(which(targetpop & !foundpop))
    tn=length(which(!targetpop & !foundpop))
    precision[i]=tp/(tp+fp)
    recall[i]=tp/(tp+fn)
    fmeasure[i]=(1+beta^2)*precision[i]*recall[i]/(beta^2*precision[i]+recall[i])
  }
  if (length(which(!is.nan(fmeasure)))==0) return(NA)
  return(gate[[which.max(fmeasure)]])
}

##Given a cell population (targetpop) and a data matrix (x) and a set of polygon gates (gates), calculate F-measures (using the given beta) of all the gates for the cells identified from the previous step (marked in currentpop).
.getScores <-
  function(x, targetpop, gates, skipdims=NA, currentpop=rep(TRUE,length(targetpop)), beta, randomize, predimx=NULL, predimy=NULL){
    if (length(which(!is.na(unlist(gates))))==0)
      return(NA)
    n=dim(x)[2]
    fmeasure=matrix(0,n,n)
    precision=matrix(0,n,n)
    recall=matrix(0,n,n)
    for (i in 1:(n-1)){
      for (j in (i+1):n){
        if (!is.null(predimx))
          if (i!=predimx)
            next
        if (!is.null(predimy))
          if (j!=predimy)
            next
        if (i %in% skipdims)
          next
        if (j %in% skipdims)
          next
        if (is.na(gates[[i]][[j]]) | (length(gates[[i]][[j]])==0)){
          fmeasure[i,j]=0
          precision[i,j]=0
          recall[i,j]=0
          next
        }
        foundpop=currentpop & as.logical(.inpoly(x[,i],x[,j],list(x=gates[[i]][[j]][,1],y=gates[[i]][[j]][,2])))
        tp=length(which(targetpop & foundpop))
        fp=length(which(!targetpop & foundpop))
        fn=length(which(targetpop & !foundpop))
        tn=length(which(!targetpop & !foundpop))
        precision[i,j]=tp/(tp+fp)
        recall[i,j]=tp/(tp+fn)
        fmeasure[i,j]=(1+beta^2)*precision[i,j]*recall[i,j]/(beta^2*precision[i,j]+recall[i,j])
      }
    }
    if (randomize){
      answer=sample(as.vector(fmeasure), size=1, prob=as.vector(fmeasure))
    }
    if (!randomize){
      answer=max(fmeasure)
    }
    temp=which(fmeasure==answer,arr.ind=TRUE)
    if((dim(temp))[1]>1)
      temp=temp[1,]
    i=temp[1]
    j=temp[2]
    newpop=rep(FALSE, length(currentpop))
    if (length(gates[[i]][[j]])!=0)
      newpop=currentpop & as.logical(.inpoly(x[,i],x[,j],list(x=gates[[i]][[j]][,1],y=gates[[i]][[j]][,2])))
    return(list(i=i,j=j,fmeasure=fmeasure[i,j],precision=precision[i,j],recall=recall[i,j],newpop=newpop,allfmeasures=fmeasure))
  }

##plot the F-measures, precisions, and recalls of the selected gating steps.
setMethod("plot", signature(x="GatingProjection"), function(x, legendpos='topright', beta=1, targetpop=NULL, lwd=2, ...){
  plot(0,col=par('bg'),xlim=c(0,length(x@fmeasure)),ylim=c(0,1),xlab='Gating Step', ylab='F-measure/Purity/Yield',col.lab=par('fg'),col.axis=par('fg'))
  if (is.null(targetpop)){
    precision=0
    fmeasure=0
    recall=1
  }
  if (!is.null(targetpop)){
    recall=1
    precision=length(which(targetpop))/length(targetpop)
    fmeasure=(1+beta^2)*precision*recall/(beta^2*precision+recall)
  }
  lines(0:length(x@fmeasure),c(fmeasure,x@fmeasure),col='red',lwd=lwd)
  lines(0:length(x@fmeasure),c(precision,x@precision),col='green',lwd=lwd)
  lines(0:length(x@fmeasure),c(recall,x@recall),col='blue',lwd=lwd)
  legend(legendpos,col=c('red','green','blue'),legend=c('F-measure','Purity', 'Yield'), lwd=lwd*1.5, inset=0.05)
})

setMethod("show", signature(object="GatingProjection"), function(object){
  cat("Object of class 'GatingProjection'","\n")
  cat("This object has the following slots: \n")
  cat("dimx, dimy, fmeasure, fmeasures, gates, pops, precision, recall, subsampleindex\n")
  cat("The selected gating strategy uses the following pairs of markers:\n")
  for (i in seq(length(object@dimx))){
    cat(sprintf('(%s, %s)',names(object@dimx)[i],names(object@dimy)[i]))
    if (i < length(object@dimx))
      cat(' -> ')
  }
  cat("\n")
})



##create a scatter plot for each of the gating steps. ncolrow is a vector of length 2 indicating the required number of rows and columns in the plot. cexs and cols: cex and color values for 1-previously excluded cells 2-non-selected cells 3-selected cells.
plot.GateFinder <- function(x, y, ncolrow=c(1,max(targetpop)), targetpop=NULL, beta=NULL,  cexs=NULL, cols=NULL, subsample=length(targetpop), max.iter=length(y@gates), pot=TRUE, xlim=NULL, ylim=NULL, asinh.axis=FALSE,...){
  if (subsample<length(targetpop)){
    subsampleindex <- sample(seq(nrow(x)),subsample)
    x=x[subsampleindex,]
    targetpop=targetpop[subsampleindex]
    for (i in 1:length(y@pops)){
      y@pops[[i]]=y@pops[[i]][subsampleindex]
    }
  }
  if (is.null(cexs))
    cexs=rep(1,3)
  if (is.null(cols)){
    cols=c('gray', par('fg'), 'red')
    if (par('bg')=='black'){
      cols=c(colors()[190], 'white', 'red')
    }
  }
  par(mfrow=ncolrow)
  for (i in 1:length(y@gates)){
    if (i>max.iter)
      next
    myxlim=xlim
    myylim=ylim
    if (is.null(xlim)){
      myxlim=c(min(x[,y@dimx[i]]),max(x[,y@dimx[i]]))
      if (asinh.axis){
        myxlim=asinh(0.2*c(-20,10000))            
      }
    }
    if (is.null(ylim)){
      myylim=c(min(x[,y@dimy[i]]),max(x[,y@dimy[i]]))
      if (asinh.axis){
        myylim=asinh(0.2*c(-20,10000))            
      }
    }
    myvect=rep(TRUE,length(targetpop))
    if (i>1)
      myvect=y@pops[[i-1]] 
    if (length(which(!myvect))==0)
      plot(0,col=par('bg'),xlim=myxlim,ylim=myylim,col.lab=par('fg'),col.axis=par('fg'), xlab=colnames(x)[y@dimx[i]], ylab=colnames(x)[y@dimy[i]], axes=FALSE)
    if (length(which(!myvect))>0)
      plot(x[!myvect,c(y@dimx[i],y@dimy[i])],pch='.',col=cols[1],xlim=myxlim,ylim=myylim,col.lab=par('fg'),col.axis=par('fg'), cex=cexs[1],axes=FALSE)
    if (pot){            
      points(x[myvect & (!targetpop),c(y@dimx[i],y@dimy[i])],pch='.',col=cols[2], cex=cexs[2])
      points(x[myvect & targetpop,c(y@dimx[i],y@dimy[i])],pch='.',col=cols[3], cex=cexs[3])
    }
    if (asinh.axis){
      ti=c(-20, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
      ti2=c(-10, -1, 1, 10, 100, 1000, 10000)
      axis(1, at=asinh(0.2*ti), labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
      axis(2, at=asinh(0.2*ti), labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
      axis(1, at=asinh(0.2*ti2), labels=ti2, las=2, col.axis=par('fg'),mgp=c(3,0.75,0), cex.axis=0.75)
      axis(2, at=asinh(0.2*ti2), labels=ti2, las=2, col.axis=par('fg'),mgp=c(3,0.75,0), cex.axis=0.75)
    }
    if (!asinh.axis){
      axis(1,col.axis=par('fg'));
      axis(2,col.axis=par('fg'));            
    }
    box()
    if (!pot){
      points(x[myvect ,c(y@dimx[i],y@dimy[i])],pch='.',col=cols[targetpop[myvect]+2], cex=cexs[targetpop[myvect]+2])
    }
    temp=y@gates[[i]]
    temp=rbind(temp,temp[1,])
    lines(temp,lwd=2,col='blue')
  }
}


##create a scatter plot for each of the gating steps. ncolrow is a vector of length 2 indicating the required number of rows and columns in the plot. cexs and cols: cex and color values for 1-previously excluded cells 2-non-selected cells 3-selected cells.
setMethod("plot", signature(x='matrix', y="GatingProjection"), function(x, y, ncolrow=c(1,max(targetpop)), targetpop=NULL, beta=NULL,  cexs=NULL, cols=NULL, subsample=length(targetpop), max.iter=length(y@gates), pot=TRUE, xlim=NULL, ylim=NULL, asinh.axis=FALSE,...){
  plot.GateFinder(x=x, y=y, ncolrow=ncolrow, targetpop=targetpop, beta=beta, cexs=cexs, cols=cols, subsample=subsample, max.iter=max.iter, pot=pot, xlim=xlim, ylim=ylim, asinh.axis=asinh.axis,...)
})



##given a set of results, apply them to a new larger dataset
.gapply <- function(x, ans, newtargetpop, beta=1){
  n=length(ans@dimx)
  pops=list()
  currentpop=rep(TRUE,dim(x)[1])
  fmeasure=vector()
  recall=vector()
  precision=vector()
  for (i in 1:length(ans@gates)){
    temp=list()
    temp$x=ans@gates[[i]][,1]
    temp$y=ans@gates[[i]][,2]
    pops[[i]]=currentpop & as.logical(.inpoly(x[,ans@dimx[i]],x[,ans@dimy[i]],temp))
    currentpop=pops[[i]]
    foundpop=currentpop
    tp=length(which(newtargetpop & foundpop))
    fp=length(which(!newtargetpop & foundpop))
    fn=length(which(newtargetpop & !foundpop))
    tn=length(which(!newtargetpop & !foundpop))
    precision[i]=tp/(tp+fp)
    recall[i]=tp/(tp+fn)
    fmeasure[i]=(1+beta^2)*precision[i]*recall[i]/(beta^2*precision[i]+recall[i])
  }
  ans2=ans
  ans2@pops=pops
  ans2@fmeasure=fmeasure
  ans2@precision=precision
  ans2@recall=recall
  return(ans2)
}

.inpoly <- function (x, y, POK) 
{
  kin = inout(cbind(x, y), cbind(POK$x, y = POK$y), 
              bound = TRUE)
  G = rep(0, length(x))
  G[kin] = 1
  return(G)
}

##rank this subspace for all target clusters
.scoreSubSpace <- function(x,targetclusters,beta=beta){
  gas=list()
  for (c in seq(max(targetclusters))){
    targetpop <- (targetclusters==c)
    ans=GateFinder(x, targetpop, max.iter=2, outlier.percentile=c(0, 0.01,0.25,0.5,0.75), subsample=round(nrow(x)*0.9))
    gas[[c]]=ans
  }
  fm=matrix(NA,max(targetclusters),max(targetclusters))
  for (i in seq(length(gas))){
    for (j in seq(length(gas))){
      foundpop=gas[[i]]@pops[[2]]
      targetpop=(targetclusters==j)
      tp=length(which(targetpop & foundpop))
      fp=length(which(!targetpop & foundpop))
      fn=length(which(targetpop & !foundpop))
      tn=length(which(!targetpop & !foundpop))
      pr=tp/(tp+fp)
      re=tp/(tp+fn)
      if (pr+re==0){
        fm[i,j]=0
        next;
      }            
      fm[i,j]=(1+beta^2)*pr*re/(beta^2*pr+re)        
    }
  }
  for (i in seq(length(gas))){
    for (j in seq(length(gas))){
      if (i>j){
        fm[j,i]=mean(c(fm[i,j],fm[j,i]))
        fm[i,j]=NA
      }
    }
  }
  fm=matrix(0,max(targetclusters),max(targetclusters))
  for (i in seq(length(gas))){
    for (j in seq(length(gas))){
      for (k in seq(length(gas))){
        foundpop=gas[[i]]@pops[[2]]
        targetpop=(targetclusters==j)
        tp=length(which(targetpop & foundpop))
        fp=length(which(!targetpop & foundpop))
        fn=length(which(targetpop & !foundpop))
        tn=length(which(!targetpop & !foundpop))
        pr=tp/(tp+fp)
        re=tp/(tp+fn)
        if (pr+re==0){
          next;
        }            
        ##fm[i,j]=(1+beta^2)*pr*re/(beta^2*pr+re)
        fm[i,j]=pr
      }
    }
  }
  up=0
  down=0
  for (i in seq(nrow(fm))){
    up=up+fm[i,i]
    down=down+sum(fm[-i])
  }
  score=up##/down
  return(score/length(gas))
}


.ah2sp <- function(x, increment=360, rnd=10, proj4string=CRS(as.character(NA))){
  if (class(x) != "ahull"){
    stop("x needs to be an ahull class object")
  }
  ## Extract the edges from the ahull object as a dataframe
  xdf <- as.data.frame(x$arcs)
  ## Remove all cases where the coordinates are all the same      
  xdf <- subset(xdf,xdf$r > 0)
  res <- NULL
  if (nrow(xdf) > 0){
    ## Convert each arc to a line segment
    linesj <- list()
    prevx<-NULL
    prevy<-NULL
    j<-1
    for(i in 1:nrow(xdf)){
      rowi <- xdf[i,]
      v <- c(rowi$v.x, rowi$v.y)
      theta <- rowi$theta
      r <- rowi$r
      cc <- c(rowi$c1, rowi$c2)
      ## Arcs need to be redefined as strings of points. Work out the number of points to allocate in this arc segment.
      ipoints <- 2 + round(increment * (rowi$theta / 2),0)
      ## Calculate coordinates from arc() description for ipoints along the arc.
      angles <- anglesArc(v, theta)
      seqang <- seq(angles[1], angles[2], length = ipoints)
      x <- round(cc[1] + r * cos(seqang),rnd)
      y <- round(cc[2] + r * sin(seqang),rnd)
      ## Check for line segments that should be joined up and combine their coordinates
      if (is.null(prevx)){
        prevx<-x
        prevy<-y
      } else if (x[1] == round(prevx[length(prevx)],rnd) && y[1] ==
                 round(prevy[length(prevy)],rnd)){
        if (i == nrow(xdf)){
          ##We have got to the end of the dataset
          prevx<-append(prevx,x[2:ipoints])
          prevy<-append(prevy,y[2:ipoints])
          prevx[length(prevx)]<-prevx[1]
          prevy[length(prevy)]<-prevy[1]
          coordsj<-cbind(prevx,prevy)
          colnames(coordsj)<-NULL
          ## Build as Line and then Lines class
          linej <- Line(coordsj)
          linesj[[j]] <- Lines(linej, ID = as.character(j))
        } else {
          prevx<-append(prevx,x[2:ipoints])
          prevy<-append(prevy,y[2:ipoints])
        }
      } else {
        ## We have got to the end of a set of lines, and there are several such sets, so convert the whole of this one to a line segment and reset.
        prevx[length(prevx)]<-prevx[1]
        prevy[length(prevy)]<-prevy[1]
        coordsj<-cbind(prevx,prevy)
        colnames(coordsj)<-NULL
        ## Build as Line and then Lines class
        linej <- Line(coordsj)
        linesj[[j]] <- Lines(linej, ID = as.character(j))
        j<-j+1
        prevx<-NULL
        prevy<-NULL
      }
    }
    ## Promote to SpatialLines
    lspl <- SpatialLines(linesj)
    ## Convert lines to polygons
    ## Pull out Lines slot and check which lines have start and end points that are the same
    lns <- slot(lspl, "lines")
    polys <- sapply(lns, function(x) { 
      crds <- slot(slot(x, "Lines")[[1]], "coords")
      identical(crds[1, ], crds[nrow(crds), ])
    }) 
    ## Select those that do and convert to SpatialPolygons
    polyssl <- lspl[polys]
    list_of_Lines <- slot(polyssl, "lines")
    sppolys <- SpatialPolygons(list(Polygons(lapply(list_of_Lines, function(x) {
      Polygon(slot(slot(x, "Lines")[[1]], "coords")) }), ID = "1")),
      proj4string=proj4string)
    ## Create a set of ids in a dataframe, then promote to
    SpatialPolygonsDataFrame
    hid <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "ID"))
    areas <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "area"))
    df <- data.frame(hid,areas)
    names(df) <- c("HID","Area")
    rownames(df) <- df$HID
    res <- SpatialPolygonsDataFrame(sppolys, data=df)
    res <- res[which(res@data$Area > 0),]
  }  
  return(res)
}

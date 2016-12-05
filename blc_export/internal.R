library(splancs)
library(igraph)

mcgaussprec <- function(pts, sds, xlim = c(0,1), ylim=c(0,1), psd=function(sd){0},minsd=0.1,maxsd=100, grid=100){
  N = dim(pts)[1]
  fsd <- Vectorize(function(sd){
    wts = 1/(sd^2 + sds^2); tildeN = sum(wts)
    mu = c(sum(wts * pts[,1])/tildeN, sum(wts*pts[,2])/tildeN)
    totdist = sum(c(wts*(pts[,1] - mu[1])^2, wts*(pts[,2]-mu[2])^2))    

    ##x-axis
    log(pnorm(sqrt(tildeN) * (xlim[2]-mu[1])) - pnorm(sqrt(tildeN) * (xlim[1] - mu[1])))+
      ##y-axis
      log(pnorm(sqrt(tildeN) * (ylim[2]-mu[2])) - pnorm(sqrt(tildeN) * (ylim[1] - mu[2])))+
        ##marginalised (with factor 2pi taken from standardisation above)
        -(N-1)*log(2 * pi)+sum(log(wts))-totdist/2+
          ##size of area
          -log(diff(xlim)*diff(ylim))+
            ##cst -- left in standardisation above
            -log(tildeN)+
              ##prior on sd
              psd(sd)
  })
  ##discrete prior:
  x = seq(minsd, maxsd, length=grid)[-1];
  values = fsd(x); dx=x[2]-x[1]; m = max(values); int=sum(exp(values-m))*dx;
  log(int)+m
}


mkcols <- function(labels){
  t=table(labels)
  cnames=names(t[t>1])  
  colors=sample(rainbow(length(cnames)))
  s=sapply(labels, function(l){
    i=which(names(t)==l);
    if (t[i]==1){"grey"}
    else {colors[which(cnames==l)]}
  })
  s  
}

toroid <- function(pts, xlim, ylim, range){
  xd=xlim[2]-xlim[1]; yd=ylim[2]-ylim[1]
  R=pts[pts[,1] >= (xlim[2]-range),,drop=FALSE]; Rshift=t(apply(R, 1, function(v) {v - c(xd, 0)}))
  L=pts[pts[,1] <= range,,drop=FALSE]; Lshift=t(apply(L, 1, function(v) {v + c(xd, 0)}))
  U=pts[pts[,2] >= ylim[2]-range,,drop=FALSE]; Ushift=t(apply(U, 1, function(v) {v - c(0, yd)}))
  D=pts[pts[,2] <= range,,drop=FALSE]; Dshift=t(apply(D, 1, function(v) {v + c(0, yd)}))

  LU=pts[(pts[,1] <= range) & (pts[,2] >= ylim[2]-range),,drop=FALSE]; LUshift=t(apply(LU, 1, function(v) {v + c(xd, -yd)}))
  RU=pts[(pts[,1] >= xlim[2]-range) & (pts[,2] >= ylim[2]-range),,drop=FALSE]; RUshift=t(apply(RU, 1, function(v) {v + c(-xd, -yd)}))
  RD=pts[(pts[,1] >= xlim[2]-range) & (pts[,2] <= range),,drop=FALSE]; RDshift=t(apply(RD, 1, function(v) {v + c(-xd, yd)}))
  LD=pts[(pts[,1] <= range) & (pts[,2] <= range),,drop=FALSE]; LDshift=t(apply(LD, 1, function(v) {v + c(xd, yd)}))
  if (length(Rshift)>0)
    pts=rbind(pts, Rshift)
  if (length(Lshift)>0)
    pts=rbind(pts, Lshift)
  if (length(Ushift)>0)
    pts=rbind(pts, Ushift)
  if (length(Dshift)>0)
    pts=rbind(pts, Dshift)
  if (length(LUshift)>0)
    pts=rbind(pts, LUshift)
  if (length(RUshift)>0)
    pts=rbind(pts, RUshift)
  if (length(RDshift)>0)
    pts=rbind(pts, RDshift)
  if (length(LDshift)>0)
    pts=rbind(pts, LDshift)
  pts
}



Kclust <- function(pts, sds=0, xlim, ylim, psd=NULL, minsd=NULL, maxsd=NULL, useplabel=TRUE, alpha=NULL, pb=.5, rseq=seq(10, 200, by=5), thseq=seq(5, 500, by=5), score=FALSE, rlabel=FALSE, report=TRUE, clustermethod="K"){
    N = dim(pts)[1]    
    if (N==1){
          rs=c()
          ths=c()
          for (r in rseq){
              for (th in thseq){
                  rs=c(rs, r); ths=c(ths,th)
              }
          }
          labels=rep(1,length(rs)); dim(labels)=c(length(rs),1)
          return(list(scores=rep(0,length(rs)), scale=rs, thresh=ths, labels=labels))
    }
    tor=toroid(pts, xlim, ylim, max(rseq))
    D=as.matrix(dist(tor))
    D=D[1:N, 1:N]
    scores=c()
    retlabels=c()
    rs=c()
    ths=c()
    for (r in rseq){
        K=apply(D, 1, function(v){sum(v <= r)-1})
        L=sqrt((diff(xlim)+2*max(rseq))*(diff(ylim)+2*max(rseq))* K /(pi * (dim(tor)[1]-1)))
        for (th in thseq){
            C=which(L>=th)
            if (clustermethod=="K"){
                if (length(C)>0){
                    G = graph.adjacency(D[C,C] < 2 *r)
                    lab=clusters(G, "weak")
                    labels=(N+1):(2*N); labels[C]=lab$membership
                }
                else labels=1:N
            }
            if (clustermethod=="DBSCAN"){
                if (length(C)>0){
                    G = graph.adjacency(D[C,C] < r)
                    lab=clusters(G, "weak")
                    labels=(N+1):(2*N); labels[C]=lab$membership
                    ##hoovering up boundary points by (arbitrarily) assigning to the first clustered
                    for (i in (1:N)[-C]){
                        closeto=which(D[C,i]< r)
                        if (length(closeto)>0) labels[i]=labels[C[closeto[1]]]
                    }
                }
                else labels=1:N            
            }            
            s=0
            if (score){
                s=scorewprec(labels=labels, pts=pts, sds=sds, xlim=xlim, ylim=ylim, psd=psd, minsd=minsd, maxsd=maxsd, useplabel=useplabel, alpha=alpha, pb=pb)
                scores=c(scores, s)
            }
            rs=c(rs, r); ths=c(ths,th); 
            if (report){cat("Scale:", r, "Thr:", th, "Score: ", s, "\n")}
            if (rlabel){
                retlabels=rbind(retlabels, labels)
            }
        }
    }
    list(scores=scores, scale=rs, thresh=ths, labels=retlabels)
}

writeRes <- function(res, rfile, labdir, bestonly=FALSE){
  scale=unique(res[["scale"]]); scale=scale[order(as.numeric(scale))]
  thresh = unique(res[["thresh"]]); thresh=thresh[order(as.numeric(thresh))]
  cat("0", scale, sep="\t", file=rfile); cat("\n", file=rfile, append=TRUE)
  for (line in thresh){
    scales=res[["scale"]][res[["thresh"]]==line]; o=order(scales); scales=scales[o]
    scores=res[["scores"]][res[["thresh"]]==line]; scores=scores[o]
    cat(line, "\t", sep="", file=rfile, append=TRUE); cat(scores, sep="\t", append=TRUE, file=rfile); cat("\n", file=rfile, append=TRUE)
  }
  dir.create(labdir)
  if (bestonly) is=which.max(res[["scores"]])
  else is=(1:dim(res[["labels"]])[1])
  for (i in is){
    f=file.path(paste(labdir, "/clusterscale", res[["scale"]][i], "\ thresh", res[["thresh"]][i], "labels.txt", sep=""))
    cat(res[["labels"]][i,], file=f, sep=","); cat("\n", file=f, append=TRUE)
  }
}

nClusters <- function(labels){
  sum(table(labels)>1)
}

percentageInCluster <- function(labels){
  Nb=sum(table(labels)==1) 
  (length(labels)-Nb)/length(labels) * 100
}

molsPerCluster <- function(labels){
  ta=table(labels); ta[ta>1]
}

nMolsPerCluster <- function(labels){
  length(labels)*percentageInCluster(labels)/(100*nClusters(labels))
}

histnMols <- function(labels){
  ta=table(labels)[table(labels)>1]; h=hist(ta, plot=FALSE)
  plot(h, xlab="Number of molecules", ylab="Number of clusters", main="")
}

clusterRadii <- function(pts, labels){
  radii=tapply(1:(dim(pts)[1]), labels, function(v){
    if (length(v)==1) -1
    else {
        mean(c(sd(pts[v,1]), sd(pts[v,2])))
    }
    })
  radii[radii>=0]
}

clusterStatistics <- function(pts, labels){
    iscluster=table(labels)>1
    if (sum(iscluster)==0) return(-1)
    clusters=names(which(iscluster))
    sapply(clusters, function(l){
        ptsl=pts[labels==l,]
        v=c(colMeans(ptsl[,1:2]), mean(c(sd(ptsl[,1]), sd(ptsl[,2]))), dim(ptsl)[1])
        dim(v)=c(4,1)
        v
    })
}

convexHullAreas <- function(pts, labels){
  areas=tapply(1:(dim(pts)[1]), labels, function(v){
    if (length(v)==1) -1
    else {
      i<-chull(pts[v,1],pts[v,2])
      areapl(as.matrix(pts[v[i],]))
    }
    })
  areas[areas>=0]
}


reldensity <- function(pts, labels, xlim, ylim){
    rs=clusterRadii(pts, labels)
    tb=table(labels)
    nclustered = sum(tb[tb>=2])
    nb=length(labels)-nclustered
    areaclustered=sum(pi*rs^2)
    (nclustered/areaclustered)/(nb/(diff(xlim)*diff(ylim)-areaclustered))
}

plabel <- function(labels, alpha, pb){
  cnt <-tapply(1:length(labels), labels, length)
  cl =cnt[cnt!=1]
  B = length(labels)-sum(cl)
  Bcont = B*log(pb)+(1-B)*log(1-pb)
  ## Green 2001 p.357, Scand J Statist 28
  partcont=0
  if (length(cl) >0)
    partcont=length(cl)*log(alpha)+lgamma(alpha)+sum(lgamma(cl))-lgamma(alpha+sum(cl))
  Bcont+partcont
}

scorewprec <- function(labels, pts, sds, xlim, ylim, psd, minsd, maxsd, useplabel=TRUE, alpha=NULL, pb=.5){  
  s=sum(tapply(1:(dim(pts)[1]), labels, function(v){
    if (length(v)>1) mcgaussprec(pts[v,], sds[v], xlim, ylim, psd=psd, minsd=minsd, maxsd=maxsd)
    else -log(diff(xlim)*diff(ylim))      
  }))  
  prlab=0
  if (useplabel){
    if (is.null(alpha)){
      cnt <-tapply(1:length(labels), labels, length)
      n =sum(cnt[cnt!=1])
      alpha=20
    }
    prlab=plabel(labels, alpha, pb)
  }
  s+prlab
}


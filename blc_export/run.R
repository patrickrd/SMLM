source("internal.R")
foldernames=c("ROIs")
sapply(foldernames, function(foldername){
r = readLines(con=file.path(paste(foldername, "/config.txt", sep="")))
get <- function(type){i = grep(type,r); strsplit(r[i], "=")[[1]][2]}
as.v <- function(ch){as.numeric(strsplit(ch,",")[[1]])}
model=get("model")
{if (model=="Gaussian(prec)"){
  xlim = as.v(get("xlim"))
  ylim = as.v(get("ylim"))
  histbins = as.v(get("histbins"))
  histvalues = as.v(get("histvalues"))
  if (length(grep("pbackground",r))==0 | length(grep("alpha",r))==0){
    useplabel=FALSE; pb=NULL; alpha=NULL
  }
  else {
    useplabel=TRUE; 
    pb=as.numeric(get("pbackground"))
    alpha=as.numeric(get("alpha"))
  }
  if (length(grep("bestonly",r))==0) bestonly=FALSE
  else bestonly=as.numeric(get("bestonly"))>0
  if (length(grep("rseq",r))==0) rseq=seq(10, 200, by=5)
  else {
      rparams=as.v(get("rseq"))
      rseq=seq(rparams[1], rparams[2], by=rparams[3])
  }
  if (length(grep("thseq",r))==0) thseq=seq(5, 500, by=5)
  else {
      thparams=as.v(get("thseq"))
      thseq=seq(thparams[1], thparams[2], by=thparams[3])
  }
  if (length(grep("clustermethod",r))==0) clustermethod="K"
  else {
      method=as.numeric(get("clustermethod"))
      if (method==1) clustermethod="K"
      else clustermethod="DBSCAN"
  }
}
else {stop("Haven't implemented anything else!")}}

o = order(histbins); histbins=histbins[o]; histvalues=histvalues[o]
f = approxfun(histbins, histvalues, yleft=histvalues[1], yright=histvalues[length(histvalues)])
cst=integrate(f, lower=histbins[o],upper=histbins[length(histbins)])$value
psd <- function(sd){
  log(f(sd))-log(cst) 
}
minsd = histbins[1]; maxsd = histbins[length(histbins)]

ld=list.dirs(foldername, recursive=FALSE)
ld=ld[ld!=foldername]

sapply(file.path(ld), function(foldername){



data= read.csv(file.path(paste(foldername, "/data.txt", sep="")))
pts = data[,1:2]; sds = data[,3];
res=Kclust(pts=pts, sds=sds, xlim=xlim, ylim=ylim, psd=psd, minsd=minsd, maxsd=maxsd, useplabel=useplabel, alpha=alpha, pb=pb, score=TRUE, rlabel=TRUE, report=TRUE, rseq=rseq, thseq=thseq, clustermethod=clustermethod)
writeRes(res, file.path(paste(foldername, "/r_vs_thresh.txt", sep="")), file.path(paste(foldername, "/labels", sep="")), bestonly=bestonly)

})


})

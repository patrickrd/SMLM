foldername="ROIs"
r = readLines(con=file.path(paste(foldername, "/sim_params.txt", sep="")))
get <- function(type){
    i = grep(type,r)
    if (length(i)==0) NA
    else strsplit(r[i], "=")[[1]][2]
}
as.v <- function(ch){as.numeric(strsplit(ch,",")[[1]])}

##multimer: if greater than 1, ignore all clustering parameters (e.g. nclusters, sdcluster) and do CSR with some proportion multimered. If zero or unspecified, normal simulation
##sdcluster: if length is one, get nclusters from nclusters, otherwise get nclusters from sdcluster
##ab: background density is fixed in any vertical line, and a Beta with parameters a and b going left to right


xlim = as.v(get("xlim"))
ylim = as.v(get("ylim"))
gammaparams=as.v(get("gammaparams"))
nsim=as.numeric(get("nsim"))

multimer=as.numeric(get("multimerisation"))

if (is.na(multimer)) multimer=0

{if (multimer <=1){
    multimer=1
    background=as.numeric(get("background"))
    sdcluster=as.v(get("sdcluster"))
    if (length(sdcluster)==1) {
        nclusters=as.numeric(get("nclusters"))
        sdcluster=rep(sdcluster, nclusters)
    }
    else nclusters=length(sdcluster)
    
    molspercluster=as.numeric(get("molspercluster"))


}
else {
    propmultimered=as.numeric(get("propmultimered"))
    nmols=as.numeric(get("nmols_for_multimer_case"))
    nclusters=floor(propmultimered*nmols)
}}

ab=get("ab")
{if (!is.na(ab)){ab=as.v(ab); fa=ab[1]; fb=ab[2]; ab=TRUE}
else ab=FALSE}


sapply(1:nsim, function(expi){
    clusteredpts <- function(n, sdcluster){
        centre=c(runif(1, min=xlim[1], max=xlim[2]), runif(1, min=ylim[1], max=ylim[2]))
        cbind(rnorm(n, mean=centre[1], sd=sdcluster), rnorm(n, mean=centre[2], sd=sdcluster))
    }
    ptsc=c()
    lc=c()
    if (multimer<=1){        
        for (i in 1:nclusters){
            lc=c(lc, rep(i,molspercluster))
            ptsc=rbind(ptsc,clusteredpts(molspercluster, sdcluster[i]))
        }
    }
    else {
        if (nclusters>0){
            for (i in 1:nclusters){
                lc=c(lc, rep(i,multimer))
                ptsc=rbind(ptsc,matrix(rep(c(runif(1, min=xlim[1], max=xlim[2]), runif(1, min=ylim[1], max=ylim[2])), each=multimer), ncol=2))
            }
        }
    }
    
    sds=rgamma(dim(ptsc)[1], gammaparams[1], gammaparams[2])
    noise=cbind(rnorm(dim(ptsc)[1]), rnorm(dim(ptsc)[1]))
    ptsn=ptsc+noise*sds
    
    inside = (ptsn[,1] >= xlim[1]) & (ptsn[,1] <= xlim[2]) & (ptsn[,2] >= ylim[1]) & (ptsn[,2] <= ylim[2])
    ptsn=ptsn[inside,]
    lc=lc[inside]
    sds=sds[inside]
    npts=dim(ptsn)[1]
    if (multimer <= 1){
        nb=ceiling(npts*background/(1-background))
    }
    else {
        nb=nmols-nclusters ##not dealing with the very unlikely case that all multimers are moved off the ROI due to localisation imprecision
    }
    
    if (!ab) ptsb=cbind(runif(nb, min=xlim[1], max=xlim[2]), runif(nb, min=ylim[1], max=ylim[2]))
    else {
        xpts=rbeta(nb, fa, fb)*diff(xlim)+xlim[1]
        ptsb=cbind(xpts[1:nb], runif(nb, min=ylim[1], max=ylim[2]))
    }
    pts=rbind(ptsn, ptsb)
    labels=1:(dim(pts)[1]); labels[1:length(lc)]=lc
    sdsb=rgamma(dim(pts)[1]-length(lc), gammaparams[1], gammaparams[2])
    sds=c(sds, sdsb)
    
    data=cbind(pts, sds, labels)
    colnames(data)=c("x","y","sd","clusterID")

    dir.create(file.path(paste(foldername, "/", expi, sep="")))
    write.csv(data, file=file.path(paste(foldername, "/", expi, "/data.txt", sep="")), row.names=FALSE, quote=FALSE)

})

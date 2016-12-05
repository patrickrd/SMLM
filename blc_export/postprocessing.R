source("internal.R")
foldernames=c("ROIs")

sapply(foldernames, function(expname){
    nexpname=expname
    r = readLines(con=file.path(paste(nexpname, "/config.txt", sep="")))
    get <- function(type){i = grep(type,r); strsplit(r[i], "=")[[1]][2]}
    as.v <- function(ch){as.numeric(strsplit(ch,",")[[1]])}
    if (length(grep("skeleton",r))==0) skeleton=FALSE else skeleton=as.logical(as.numeric(get("skeleton")))
    
    if (skeleton){
        nexpname=paste("R_", expname, sep="")
        dir.create(file.path(nexpname))
        file.copy(file.path(paste(expname, "/config.txt", sep="")), paste(nexpname, "/config.txt", sep=""))
        file.copy(file.path(paste(expname, "/sim_params.txt", sep="")), paste(nexpname, "/sim_params.txt", sep=""))
    }

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
    }
     else {stop("Haven't implemented anything else!")}}
    if (length(grep("makeplot",r))==0) makeplot=FALSE else makeplot=as.logical(as.numeric(get("makeplot")))
    if (length(grep("superplot",r))==0) superplot=FALSE else superplot=as.logical(as.numeric(get("superplot")))

    all=list.files(expname)
    dirnames=all[file.info(file.path(paste(expname, "/", all, sep="")))$isdir]
    axes=FALSE;
    cex=1/(3*sqrt(length(dirnames)))
    if (makeplot & superplot) {
        ##These settings control image size and resolution
        png(file.path(paste(nexpname, "/together.png",sep="")), width=10, height=10, units="cm", res=1200)
        nrow=ceiling(sqrt(length(dirnames)))
        par(mfrow=c(nrow, nrow))
    }

    res=lapply(dirnames, function(dirname){
        foldername=file.path(paste(expname, "/", dirname, sep=""))
        nfoldername=file.path(paste(nexpname, "/", dirname, sep=""))
        if (skeleton){
            dir.create(nfoldername)  
            file.copy(file.path(paste(foldername, "/data.txt", sep="")), file.path(paste(nfoldername, "/data.txt", sep="")))
        }
        data=read.csv(file.path(paste(nfoldername, "/data.txt", sep="")))
        
        pts = data[,1:2]; sds = data[,3];
        if (skeleton){
            file.copy(file.path(paste(foldername, "/r_vs_thresh.txt", sep="")), file.path(paste(nfoldername, "/r_vs_thresh.txt", sep="")))
        }
        r = read.csv(file.path(paste(nfoldername, "/r_vs_thresh.txt",sep="")), header=FALSE, sep="\t")
        
        m = as.matrix(r)
        cs=(m[1,])[-1]
        thr=(m[,1])[-1]
        m = as.matrix(m[2:length(m[,1]),2:length(m[1,])])
        which.maxm <- function(mat){
            indcol <- rep(1:ncol(mat), each=nrow(mat))[which.max(mat)] 
            indrow <- rep(1:nrow(mat), ncol(mat))[which.max(mat)]
            c(indrow, indcol)
        }
        best=which.maxm(m)
        bestcs=cs[best[2]]
        bestthr=thr[best[1]]
        bfile=file.path(paste(foldername, "/labels/clusterscale", bestcs, " thresh", bestthr, "labels.txt", sep=""))
        nbfile=bfile
        if (skeleton){
            dir.create(paste(nfoldername, "/labels", sep=""))    
            nbfile=file.path(paste(nfoldername, "/labels/clusterscale", bestcs, " thresh", bestthr, "labels.txt", sep=""))
            file.copy(bfile, nbfile)
        }
        labelsbest = strsplit(readLines(nbfile),",")[[1]]

        
        ##Some summaries
        wfile=file.path(paste(nfoldername, "/summary.txt", sep=""))
        cat("The best: clusterscale", bestcs, " thresh", bestthr, "labels.txt\nNumber of clusters:", nClusters(labelsbest), "\nPercentage in clusters: ", percentageInCluster(labelsbest), "%\nMean number of molecules per cluster: ", nMolsPerCluster(labelsbest), "\nMean radius: ", mean(clusterRadii(pts, labelsbest)), sep="", file=wfile)
        s=clusterStatistics(pts, labelsbest)
        if (!is.null(s)){
            wfile=file.path(paste(nfoldername, "/cluster-statistics.txt", sep=""))
            cat("x,y,sd,nmol\n", file=wfile)
            for (i in 1:dim(s)[2]){
                cat(s[,i], sep=",", append=TRUE, file=wfile); cat("\n",append=TRUE,file=wfile)
            }
        }        

        if (makeplot){
            if (!superplot){
                    pdf(file.path(paste(nfoldername, "/plot.pdf", sep="")))
                    axes=TRUE
                    cex=1
                }           
            if ("clusterID" %in% colnames(data) & !superplot){
                labelstrue=sapply(as.numeric(data[,4]), function(n){if (n==0) paste(runif(1)) else {paste(n)}})

                par(pty="s")
                par(mfrow=c(1,2))
                par(mar=c(4,4,.5, .5))
                par(oma=c(1,1,1,1))
                plot(pts, xlim=xlim, ylim=ylim, col=mkcols(labelstrue), sub="True labels", xlab="",ylab="")
                par(mar=c(4,4,.5, .5))
                plot(pts, xlim=xlim, ylim=ylim, col=mkcols(labelsbest), sub="Estimated",xlab="",ylab="")
                if (!superplot) dev.off()
            }
            else {
                par(pty="s")
                par(mar=c(0,0,0,0))
                plot(pts, xlim=xlim, ylim=ylim, col=mkcols(labelsbest), sub="Clustering",xlab="",ylab="", pch=16, cex=cex, axes=axes)
                box()
                if (!superplot) dev.off()
            }
        }
        list(radii=clusterRadii(pts, labelsbest), nmols=molsPerCluster(labelsbest), nclusters=nClusters(labelsbest), pclustered=percentageInCluster(labelsbest), totalmols=length(labelsbest), reldensity=reldensity(pts, labelsbest, xlim, ylim))
    })

    if (makeplot & superplot) dev.off()
    nmols=c()
    for (i in 1:length(res)){
        nmols=c(nmols, res[[i]]$nmols)
    }

    h=hist(nmols, plot=FALSE)
    pdf(file.path(paste(nexpname, "/nmols.pdf", sep="")))
    plot(h, xlab="Number of molecules", ylab="Number of clusters", main="")
    dev.off()
    f=file.path(paste(nexpname, "/nmols.txt", sep="")); cat(nmols, file=f, sep=","); cat("\n", file=f, append=TRUE)
    
    radii=c()
    for (i in 1:length(res)){
        radii=c(radii, res[[i]]$radii)
    }

    h=hist(radii, plot=FALSE)
    pdf(file.path(paste(nexpname, "/radii.pdf", sep="")))
    plot(h, xlab="Cluster radius", ylab="Number of clusters", main="")
    dev.off()
    f=file.path(paste(nexpname, "/radii.txt", sep="")); cat(radii, file=f, sep=","); cat("\n", file=f, append=TRUE)


    nclusters=c()
    for (i in 1:length(res)){
        nclusters=c(nclusters, res[[i]]$nclusters)
    }
    
    h=hist(nclusters, plot=FALSE)
    pdf(file.path(paste(nexpname, "/nclusters.pdf", sep="")))
    plot(h, xlab="Number of clusters", ylab="Number of regions", main="")
    dev.off()
    f=file.path(paste(nexpname, "/nclusters.txt", sep="")); cat(nclusters, file=f, sep=","); cat("\n", file=f, append=TRUE)

    pclustered=c()
    for (i in 1:length(res)){
        pclustered=c(pclustered, res[[i]]$pclustered)
    }
    
    h=hist(pclustered, plot=FALSE)
    pdf(file.path(paste(nexpname, "/pclustered.pdf", sep="")))
    plot(h, xlab="Percentage clustered", ylab="Number of regions", main="")
    dev.off()
    f=file.path(paste(nexpname, "/pclustered.txt", sep="")); cat(pclustered, file=f, sep=","); cat("\n", file=f, append=TRUE)   

    totalmols=c()
    for (i in 1:length(res)){
        totalmols=c(totalmols, res[[i]]$totalmols)
    }
    
    h=hist(totalmols, plot=FALSE)
    pdf(file.path(paste(nexpname, "/totalmols.pdf", sep="")))
    plot(h, xlab="Total Mols per ROI", ylab="Number of regions", main="")
    dev.off()
    f=file.path(paste(nexpname, "/totalmols.txt", sep="")); cat(totalmols, file=f, sep=","); cat("\n", file=f, append=TRUE)   


    reldensity=c()
    for (i in 1:length(res)){
        reldensity=c(reldensity, res[[i]]$reldensity)
    }
    
    h=hist(reldensity, plot=FALSE)
    pdf(file.path(paste(nexpname, "/reldensity.pdf", sep="")))
    plot(h, xlab="Total Mols per ROI", ylab="Number of regions", main="")
    dev.off()
    f=file.path(paste(nexpname, "/reldensity.txt", sep="")); cat(reldensity, file=f, sep=","); cat("\n", file=f, append=TRUE)

})

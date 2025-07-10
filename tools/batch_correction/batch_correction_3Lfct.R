# Author: jfmartin
# Modified by : mpetera
###############################################################################
# Correction of analytical effects inter and intra batch on intensities using quality control pooled samples (QC-pools)
# according to the algorithm mentioned by Van der Kloet (J Prot Res 2009).
# Parameters : a dataframe of Ions intensities and an other of samples? metadata which must contains at least the three following columns :
#   "batch" to identify the batches of analyses ; need at least 3 QC-pools for linear adjustment and 8 for lo(w)ess adjustment
#   "injectionOrder" integer defining the injection order of all samples : QC-pools and analysed samples
#   "sampleType" indicates if defining a sample with "sample" or a QC-pool with "pool"
# NO MISSING DATA are allowed
# Version 0.91 insertion of ok_norm function to assess correction feasibility 
# Version 0.92 insertion of slope test in ok_norm
# Version 0.93 name of log file define as a parameter of the correction function
# Version 0.94 Within a batch, test if all QCpools or samples values = 0. Definition of an error code in ok_norm function (see function for details)
# Version 0.99 include non linear lowess correction. 
# Version 1.00 the corrected result matrix is return transposed in Galaxy
# Version 1.01 standard deviation=0 instead of sum of value=0 is used to assess constant data in ok_norm function. Negative values in corrected matrix are converted to 0. 
# Version 1.02 plotsituation create a result file with the error code of non execution of correction set by function ok_norm
# Version 1.03 fix bug in plot with "reg" option. suppression of ok_norm=4 condition if ok_norm function
# Version 2.00 Addition of loess function, correction indicator, plots ; modification of returned objects' format, some plots' displays and ok_norm ifelse format
# Version 2.01 Correction for pools negative values earlier in norm_QCpool
# Version 2.10 Script refreshing ; vocabulary adjustment ; span in parameters for lo(w)ess regression ; conditionning for third line ACP display ; order in loess display
# Version 2.11 ok1 and ok2 permutation (ok_norm) ; conditional display of regression (plotsituation) ; grouping of linked lignes + conditioning (normX) ; conditioning for CVplot
# Version 2.20 acplight function added from previous toolBox.R [# Version 1.01 "NA"-coding possibility added in acplight function]
# Version 2.30 addition of suppressWarnings() for known and controlled warnings ; suppression of one useless "cat" message ; change in Rdata names ; 'batch(es)' in cat
# Version 2.90 change in handling of generated negative and Inf values
# Version 2.91 Plot improvement
# Version 3.00 - handling of sample tags' parameters
#              - accepting sample types beyond "pool" and "sample"
#              - dealing with NA
#              - changes in the normalisation strategy regarding mean values to adjust for NA or 0 values
#              - changes in the normalisation strategy regarding unconsistant values (negative or Inf)

ok_norm=function(qcp,qci,spl,spi,method,normref=NA,valimp="0") {
  # Function used for one ion within one batch to determine whether or not batch correction is possible
  # ok_norm values :
  #   0 : no preliminary-condition problem
  #   1 : standard deviation of QC-pools or samples = 0
  #   2 : insufficient number of QC-pools within a batch (n=3 for linear, n=8 for lowess or loess)
  #   2.5 : less than 2 samples within a batch
  #   3 : significant difference between QC-pools' and samples' means
  #   4 : denominator =0 when on 1 pool per batch <> 0
  #   5 : (linear regression only) the slopes ratio ?QC-pools/samples? is lower than -0.2
  #   6 : (linear regression only) none of the pool or sample could be corrected if negative and infinite values are turned into NA
  # Parameters:
  # qcp: intensity of a given ion for pools
  # qci: injection numbers for pools
  # spl: intensity of a given ion for samples
  # spi: injection numbers for samples
  # method: to provide specific checks for "linear"
  
  ok=0
  if (method=="linear") {minQC=3} else {minQC=8}
  if (length(qcp[!is.na(qcp)])<minQC) { ok=2 } else { if (length(spl[!is.na(spl)])<2) { ok=2.5
  } else { 
    if (sd(qcp,na.rm=TRUE)==0 | sd(spl,na.rm=TRUE)==0) { ok=1 
    } else {
    cvp= sd(qcp,na.rm=TRUE)/mean(qcp,na.rm=TRUE); cvs=sd(spl,na.rm=TRUE)/mean(spl,na.rm=TRUE)
    rttest=t.test(qcp,y=spl)
    reslsfit=lsfit(qci, qcp)
    reslsfitSample=lsfit(spl, spi)
    ordori=reslsfit$coefficients[1]
    penteB=reslsfit$coefficients[2]
    penteS=reslsfitSample$coefficients[2]
    # Significant difference between samples and pools
    if (rttest$p.value < 0.01) { ok=3 
    } else {
      # to avoid denominator =0 when on 1 pool per batch <> 0
      if (method=="linear" & length(which(((penteB*qci)+ordori)==0))>0 ){ ok=6 
      } else {
        # different sloop between samples and pools
        if (method=="linear" & penteB/penteS < -0.20) { ok=5 
        } else {
          # 
          if (method=="linear" & !is.na(normref) & valimp=="NA") {
            denom = (penteB * c(spi,qci) + ordori)
            normval = c(spl,qcp)*normref / denom
            if(length(which((normval==Inf)|(denom<1)))==length(normval)){ok=6}
          } 
  }}}}}} 
  ok_norm=ok
}

plotsituation <- function (x, nbid,outfic="plot_regression.pdf", outres="PreNormSummary.txt",fact="batch",span="none",
                           sm_meta=list(batch="batch", injectionOrder="injectionOrder", sampleType="sampleType", 
                           sampleTag=list(pool="pool",blank="blank",sample="sample"))) {
    # Checks for all ions in every batch if linear or lo(w)ess correction is possible.
    # Uses ok_norm function and creates a file (PreNormSummary.txt) with the corresponding error codes.
    # Also creates a pdf file with plots of linear and lo(w)ess regression lines.
    # Parameters:
    # x: dataframe with ions in columns and samples in rows ; x is the result of concatenation of sample metadata file and ions file
    # nbid: number of samples description columns (id and factors) with at least : "batch","injectionOrder","sampleType"
    # outfic: name of regression plots pdf file
    # outres: name of summary table file
    # fact: factor to be used as categorical variable for plots and PCA
    # span: span value for lo(w)ess regression; "none" for linear or default values
    # sm_meta: list of information about sample metadata coding
    indfact=which(dimnames(x)[[2]]==fact)
    indtypsamp=which(dimnames(x)[[2]]==sm_meta$sampleType) 
    indbatch=which(dimnames(x)[[2]]==sm_meta$batch)
    indinject=which(dimnames(x)[[2]]==sm_meta$injectionOrder)
    lastIon=dim(x)[2]
    nbi=lastIon-nbid # Number of ions = total number of columns - number of identifying columns
    nbb=length(levels(x[[sm_meta$batch]])) # Number of batch = number of levels of "batch" comlumn (factor)
    nbs=length(x[[sm_meta$sampleType]][x[[sm_meta$sampleType]] %in% sm_meta$sampleTag$sample])# Number of samples = number of rows with "sample" value in sampleType
    pdf(outfic,width=27,height=7*ceiling((nbb+2)/3))
    cat(nbi," ions ",nbb," batch(es) \n")
    cv=data.frame(matrix(0,nrow=nbi,ncol=2))# initialisation de la dataset qui contiendra les CV
    pre_bilan=matrix(0,nrow=nbi,ncol=3*nbb) # dataset of ok_norm function results
    for (p in 1:nbi) {# for each ion
        par (mfrow=c(ceiling((nbb+2)/3),3),ask=F,cex=1.2)
        labion=dimnames(x)[[2]][p+nbid]
        indpool=which(x[[sm_meta$sampleType]] %in% sm_meta$sampleTag$pool) # QCpools subscripts in x 
        pools1=x[indpool,p+nbid]; cv[p,1]=sd(pools1,na.rm=TRUE)/mean(pools1,na.rm=TRUE)# CV before correction
        for (b in 1:nbb) {# for each batch...
            xb=data.frame(x[(x[[sm_meta$batch]]==levels(x[[sm_meta$batch]])[b]),c(indtypsamp,indinject,p+nbid)])
            indpb = which(xb[[sm_meta$sampleType]] %in% sm_meta$sampleTag$pool)# QCpools subscripts in the current batch
            indsp = which(xb[[sm_meta$sampleType]] %in% sm_meta$sampleTag$sample)# samples subscripts in the current batch
            normLinearTest=ok_norm(xb[indpb,3],xb[indpb,2], xb[indsp,3],xb[indsp,2],method="linear",normref=mean(xb[c(indpb,indsp),3],na.rm=TRUE),valimp="NA")
            normLoessTest=ok_norm(xb[indpb,3],xb[indpb,2], xb[indsp,3],xb[indsp,2],method="loess")
            normLowessTest=ok_norm(xb[indpb,3],xb[indpb,2], xb[indsp,3],xb[indsp,2],method="lowess")
            pre_bilan[ p,3*b-2]=normLinearTest
            pre_bilan[ p,3*b-1]=normLoessTest
            pre_bilan[ p,3*b]=normLowessTest
          if(length(indpb)>1){ 
            if(span=="none"){span1<-1 ; span2<-2*length(indpool)/nbs}else{span1<-span ; span2<-span}
            if(normLoessTest!=2){resloess=loess(xb[indpb,3]~xb[indpb,2],span=span1,degree=2,family="gaussian",iterations=4,surface="direct")}
            if(length(which(!(is.na(xb[indsp,3]))))>1){resloessSample=loess(xb[indsp,3]~xb[indsp,2],span=2*length(indpool)/nbs,degree=2,family="gaussian",iterations=4,surface="direct") }
            if(normLowessTest!=2){reslowess=lowess(xb[indpb,2],xb[indpb,3],f=span2)}
            if(length(which(!(is.na(xb[indsp,3]))))>1){reslowessSample=lowess(xb[indsp,2],xb[indsp,3])}
            liminf=min(xb[,3],na.rm=TRUE);limsup=max(xb[,3],na.rm=TRUE)
            firstinj=min(xb[,2],na.rm=TRUE);lastinj=max(xb[,2],na.rm=TRUE)
            plot(xb[indsp,2],xb[indsp,3],pch=16, main=paste(labion,"batch ",b),ylab="intensity",xlab="injection order",ylim=c(liminf,limsup),xlim=c(firstinj,lastinj))
            if(nrow(xb)>(length(indpb)+length(indsp))){points(xb[-c(indpb,indsp),2], xb[-c(indpb,indsp),3],pch=18,col="grey")}
            points(xb[indpb,2], xb[indpb,3],pch=5)
            if(normLoessTest!=2){points(cbind(resloess$x,resloess$fitted)[order(resloess$x),],type="l",col="green3")}
            if(length(which(!(is.na(xb[indsp,3]))))>1){points(cbind(resloessSample$x,resloessSample$fitted)[order(resloessSample$x),],type="l",col="green3",lty=2)}
            if(normLowessTest!=2){points(reslowess,type="l",col="red")}; if(length(which(!(is.na(xb[indsp,3]))))>1){points(reslowessSample,type="l",col="red",lty=2)}
            abline(lsfit(xb[indpb,2],xb[indpb,3]),col="blue")
            if(length(which(!(is.na(xb[indsp,3]))))>1){abline(lsfit(xb[indsp,2],xb[indsp,3]),lty=2,col="blue")}
            legend("topleft",c("pools","samples"),lty=c(1,2),bty="n")
            legend("topright",c("linear","lowess","loess"),lty=1,col=c("blue","red","green3"),bty="n")
          } else {
            plot.new()
            legend("center","Plot only available when the\nbatch contains at least 2 pools.")
          }
        }
        # series de plot avant correction
        minval=min(x[p+nbid],na.rm=TRUE);maxval=max(x[p+nbid],na.rm=TRUE)
        plot( x[[sm_meta$injectionOrder]], x[,p+nbid],col=x[[sm_meta$batch]],ylim=c(minval,maxval),ylab=labion,
              main=paste0("before correction (CV for pools = ",round(cv[p,1],2),")"),xlab="injection order")
        suppressWarnings(plot.design( x[c(indtypsamp,indbatch,indfact,p+nbid)],main="factors effect before correction"))
    }
    dev.off()
    pre_bilan=data.frame(pre_bilan)
    labion=dimnames(x)[[2]][nbid+1:nbi]
    for (i in 1:nbb) {
      dimnames(pre_bilan)[[2]][3*i-2]=paste("batch",i,"linear")
      dimnames(pre_bilan)[[2]][3*i-1]=paste("batch",i,"loess")
      dimnames(pre_bilan)[[2]][3*i]=paste("batch",i,"lowess")
    }
    bilan=data.frame(labion,pre_bilan)
    write.table(bilan,file=outres,sep="\t",row.names=F,quote=F)
}


normlowess=function (xb,detail="no",vref=1,b,span=NULL,valneg="none",sm_meta=list(batch="batch", injectionOrder="injectionOrder", sampleType="sampleType", 
                     sampleTag=list(pool="pool",blank="blank",sample="sample")),min_norm=1){
  # Correction function applied to 1 ion in 1 batch. 
  # Uses a lowess regression computed on QC-pools in order to correct samples intensity values
  # xb: dataframe for 1 ion in columns and samples in rows.
  # vref: reference value (average of ion)
  # b: batch subscript
  # detail: level of detail in the outlog file
  # span: span value for lo(w)ess regression; NULL for default values
  # valneg: to determine what to do with generated negative and Inf values
  # sm_meta: list of information about sample metadata coding
  # min_norm: minimum value accepted for normalisation term (denominator); should be strictly positive
  indpb = which(xb[[sm_meta$sampleType]] %in% sm_meta$sampleTag$pool) # pools subscripts of current batch
  indsp = which(xb[[sm_meta$sampleType]] %in% sm_meta$sampleTag$sample) # samples of current batch subscripts
  labion=dimnames(xb)[[2]][3]
  newval=xb[[3]] # initialisation of corrected values = intial values
  ind <- 0 # initialisation of correction indicator
  normTodo=ok_norm(xb[indpb,3],xb[indpb,2], xb[indsp,3],xb[indsp,2],method="lowess")
  #cat("batch:",b," dim xb=",dim(xb)," ok=",normTodo,"\n")
  if (normTodo==0) {
    if(length(span)==0){span2<-2*length(indpb)/length(indsp)}else{span2<-span}
    reslowess=lowess(xb[indpb,2],xb[indpb,3],f=span2) # lowess regression with QC-pools 
    if(length(which(reslowess$y<min_norm))!=0){ # to handle cases where 0<denominator<min_norm or negative
      toajust <- which(reslowess$y<min_norm)
      if(valneg=="NA"){ reslowess$y[toajust] <- NA
      } else { if(valneg=="0"){ reslowess$y[toajust] <- -1
        } else {
          mindenom <- min(reslowess$y[reslowess$y>=min_norm],na.rm=TRUE)
          reslowess$y[toajust] <- mindenom
    } } }
    for(j in 1:nrow(xb)) {
      if (j %in% indpb) {
        newval[j]=(vref*xb[j,3]) / (reslowess$y[which(indpb==j)])
      } else { # for samples other than pools, the correction value "corv" correspond to the nearest QCpools 
        corv= reslowess$y[which(abs(reslowess$x-xb[j,2])==min(abs(reslowess$x-xb[j,2]),na.rm=TRUE))]
        if (length(corv)>1) {corv=corv[1]}
        newval[j]=(vref*xb[j,3]) / corv
      }
      if((!is.na(newval[j]))&(newval[j]<0)){newval[j]<-0}
    }
    if (detail=="reg") {
      liminf=min(xb[,3],na.rm=TRUE);limsup=max(xb[,3],na.rm=TRUE)
      firstinj=min(xb[,2],na.rm=TRUE);lastinj=max(xb[,2],na.rm=TRUE)
      plot(xb[indsp,2],xb[indsp,3],pch=16,main=paste(labion,"batch ",b),ylab="intensity",xlab="injection order",ylim=c(liminf,limsup),xlim=c(firstinj,lastinj))
      if(nrow(xb)>(length(indpb)+length(indsp))){points(xb[-c(indpb,indsp),2], xb[-c(indpb,indsp),3],pch=18)}
      points(xb[indpb,2], xb[indpb,3],pch=5)
      points(reslowess,type="l",col="red")
    }
    ind <- 1
  } else {# if ok_norm != 0 , we perform a correction based on batch pool or sample average
    if((length(which(!is.na(xb[indpb,3])))>0)&(length(which(xb[indpb,3]>0))>0)){
      moypool=mean(xb[indpb,3],na.rm=TRUE)
      newval = (vref*xb[,3])/moypool
    } else { 
      moysamp=mean(xb[indsp,3],na.rm=TRUE)
      if((!is.na(moysamp))&(moysamp>0)){
        cat("Warning: no pool value >0 detected in batch",b,"of ion",labion,": sample mean used as normalisation term.\n")
        newval = (vref*xb[,3])/moysamp
      } else {
        dev.off()
        stop(paste("\n- - - -\nNo pool nor sample value >0 in batch",b,"of ion",labion,"- correction process aborted.\n- - - -\n"))
      }
    }
  }
  newval <- list(norm.ion=newval,norm.ind=ind)
  return(newval)
}

normlinear <- function (xb,detail="no",vref=1,b,valneg="none",sm_meta=list(batch="batch", injectionOrder="injectionOrder", sampleType="sampleType", 
                        sampleTag=list(pool="pool",blank="blank",sample="sample")),min_norm=1){
  # Correction function applied to 1 ion in 1 batch.
  # Uses a linear regression computed on QC-pools in order to correct samples intensity values
  # xb: dataframe with ions in columns and samples in rows; x is a result of concatenation of sample metadata file and ion file 
  # detail: level of detail in the outlog file
  # vref: reference value (average of ion)
  # b: which batch it is
  # valneg: to determine what to do with generated negative and Inf values
  # sm_meta: list of information about sample metadata coding
  # min_norm: minimum value accepted for normalisation term (denominator); should be strictly positive
  indpb = which(xb[[sm_meta$sampleType]] %in% sm_meta$sampleTag$pool)# pools subscripts of current batch
  indsp = which(xb[[sm_meta$sampleType]] %in% sm_meta$sampleTag$sample)# samples of current batch subscripts
  labion=dimnames(xb)[[2]][3]
  newval=xb[[3]] # initialisation of corrected values = intial values
  ind <- 0 # initialisation of correction indicator
  normTodo=ok_norm(xb[indpb,3],xb[indpb,2], xb[indsp,3],xb[indsp,2],method="linear",normref=vref,valimp=valneg)
  if (normTodo==0) {
    ind <- 1
    reslsfit=lsfit(xb[indpb,2],xb[indpb,3])       # linear regression for QCpools 
    reslsfitSample=lsfit(xb[indsp,2],xb[indsp,3]) # linear regression for samples
    ordori=reslsfit$coefficients[1]
    pente=reslsfit$coefficients[2]
    if (detail=="reg") {
      liminf=min(xb[,3],na.rm=TRUE);limsup=max(xb[,3],na.rm=TRUE)
      firstinj=min(xb[,2],na.rm=TRUE);lastinj=max(xb[,2],na.rm=TRUE)
      plot(xb[indsp,2],xb[indsp,3],pch=16,
           main=paste(labion,"batch ",b),ylab="intensity",xlab="injection order",ylim=c(liminf,limsup),xlim=c(firstinj,lastinj))
      if(nrow(xb)>(length(indpb)+length(indsp))){points(xb[-c(indpb,indsp),2], xb[-c(indpb,indsp),3],pch=18)}
      points(xb[indpb,2], xb[indpb,3],pch=5)
      abline(reslsfit)
      abline(reslsfitSample,lty=2)
    }
    # correction with rescaling of ion global intensity (vref)
    newval = (vref*xb[,3]) / (pente * (xb[,2]) + ordori)
    newval[which((pente * (xb[,2]) + ordori)<min_norm)] <- -1 # to handle cases where 0<denominator<1 or negative
    # handling if any negative values
    if(length(which((newval==Inf)|(newval<0)))!=0){
      toajust <- which((newval==Inf)|(newval<0))
      if(valneg=="NA"){ newval[toajust] <- NA
      } else { if(valneg=="0"){ newval[toajust] <- 0
        } else {
          mindenom <- (pente * (xb[,2]) + ordori)
          mindenom <- min(mindenom[mindenom>=min_norm],na.rm=TRUE)
          newval[toajust] <- vref * (xb[,3][toajust]) / mindenom
        }
      }
    }
  } else {# if ok_norm != 0 , we perform a correction based on batch pool or sample average
    if((length(which(!is.na(xb[indpb,3])))>0)&(length(which(xb[indpb,3]>0))>0)){
      moypool=mean(xb[indpb,3],na.rm=TRUE)
      newval = (vref*xb[,3])/moypool
    } else { 
      moysamp=mean(xb[indsp,3],na.rm=TRUE)
      if((!is.na(moysamp))&(moysamp>0)){
        cat("Warning: no pool value >0 detected in batch",b,"of ion",labion,": sample mean used as normalisation term.\n")
        newval = (vref*xb[,3])/moysamp
      } else {
        dev.off()
        stop(paste("\n- - - -\nNo pool nor sample value >0 in batch",b,"of ion",labion,"- correction process aborted.\n- - - -\n"))
      }
    }
  }
  newval <- list(norm.ion=newval,norm.ind=ind)
  return(newval)
}


normloess <- function (xb,detail="no",vref=1,b,span=NULL,valneg="none",sm_meta=list(batch="batch", injectionOrder="injectionOrder", sampleType="sampleType", 
                       sampleTag=list(pool="pool",blank="blank",sample="sample")),min_norm=1){
    # Correction function applied to 1 ion in 1 batch. 
    # Uses a loess regression computed on QC-pools in order to correct samples intensity values.
    # xb: dataframe for 1 ion in columns and samples in rows.
    # detail: level of detail in the outlog file.
    # vref: reference value (average of ion)
    # b: batch subscript
    # span: span value for lo(w)ess regression; NULL for default values
    # valneg: to determine what to do with generated negative and Inf values
    # sm_meta: list of information about sample metadata coding
    # min_norm: minimum value accepted for normalisation term (denominator); should be strictly positive
    indpb = which(xb[[sm_meta$sampleType]] %in% sm_meta$sampleTag$pool) # pools subscripts of current batch
    indsp = which(xb[[sm_meta$sampleType]] %in% sm_meta$sampleTag$sample) # samples of current batch subscripts
    indbt = which(xb[[sm_meta$sampleType]] %in% c(sm_meta$sampleTag$sample,sm_meta$sampleTag$pool))# batch subscripts of samples and QCpools
    labion=dimnames(xb)[[2]][3]
    newval=xb[[3]] # initialisation of corrected values = intial values
    ind <- 0 # initialisation of correction indicator
    normTodo=ok_norm(xb[indpb,3],xb[indpb,2], xb[indsp,3],xb[indsp,2],method="loess")
    if (normTodo==0) { 
        if(length(span)==0){span1<-1}else{span1<-span}
        resloess=loess(xb[indpb,3]~xb[indpb,2],span=span1,degree=2,family="gaussian",iterations=4,surface="direct") # loess regression with QCpools 
        corv=predict(resloess,newdata=xb[,2])
        if(length(which(corv<min_norm))!=0){ # unconsistant values handling
          toajust <- which(corv<min_norm)
          if(valneg=="NA"){ corv[toajust] <- NA
          } else { if(valneg=="0"){ corv[toajust] <- -1
            } else {
              mindenom <- min(corv[corv>=min_norm],na.rm=TRUE)
              corv[toajust] <- mindenom
            }
          }
        } 
        newvalps=(vref*xb[indbt,3]) / corv[indbt] # to check if correction generates outlier values
        refthresh=max(c(3*(quantile(newvalps,na.rm=TRUE)[4]),1.3*(xb[indbt,3])),na.rm=TRUE)
        if(length(which(newvalps>refthresh))>0){ # if outliers
          # in this case no modification of initial value
          newval <- xb[,3]
        } else { 
          newval=(vref*xb[,3]) / corv 
          newval[newval<0] <- 0
          ind <- 1 # confirmation of correction 
        }
        if ((detail=="reg")&(ind==1)) { # plot
            liminf=min(xb[,3],na.rm=TRUE);limsup=max(xb[,3],na.rm=TRUE)
            firstinj=min(xb[,2],na.rm=TRUE);lastinj=max(xb[,2],na.rm=TRUE)
            plot(xb[indsp,2],xb[indsp,3],pch=16,main=paste(labion,"batch ",b),ylab="intensity",xlab="injection order",ylim=c(liminf,limsup),xlim=c(firstinj,lastinj))
            if(nrow(xb)>(length(indpb)+length(indsp))){points(xb[-c(indpb,indsp),2], xb[-c(indpb,indsp),3],pch=18)}
            points(xb[indpb,2], xb[indpb,3],pch=5)
            points(cbind(resloess$x,resloess$fitted)[order(resloess$x),],type="l",col="red")
        }
    } 
    if (ind==0) {# if ok_norm != 0 or if correction creates outliers, we perform a correction based on batch pool or sample average
      if((length(which(!is.na(xb[indpb,3])))>0)&(length(which(xb[indpb,3]>0))>0)){
        moypool=mean(xb[indpb,3],na.rm=TRUE)
        newval = (vref*xb[,3])/moypool
      } else { 
        moysamp=mean(xb[indsp,3],na.rm=TRUE)
        if((!is.na(moysamp))&(moysamp>0)){
          cat("Warning: no pool value >0 detected in batch",b,"of ion",labion,": sample mean used as normalisation term.\n")
          newval = (vref*xb[,3])/moysamp
        } else {
          dev.off()
          stop(paste("\n- - - -\nNo pool nor sample value >0 in batch",b,"of ion",labion,"- correction process aborted.\n- - - -\n"))
        }
      }
    }
    newval <- list(norm.ion=newval,norm.ind=ind)
    return(newval)
}



norm_QCpool <- function (x, nbid, outlog, fact, metaion, detail="no", NormMoyPool=FALSE, NormInt=FALSE, method="linear",span="none",valNull="0",
                         sm_meta=list(batch="batch", injectionOrder="injectionOrder", sampleType="sampleType", 
                         sampleTag=list(pool="pool",blank="blank",sample="sample")),min_norm=1) {
  ### Correction applying linear or lo(w)ess correction function on all ions for every batch of a dataframe.
  # x: dataframe with ions in column and samples' metadata
  # nbid: number of sample description columns (id and factors) with at least "batch", "injectionOrder", "sampleType"
  # outlog: name of regression plots and PCA pdf file
  # fact: factor to be used as categorical variable for plots  
  # metaion: dataframe of ions' metadata
  # detail: level of detail in the outlog file. detail="no" ACP + boxplot of CV before and after correction.
  #          detail="plot" with plot for all batch before and after correction.
  #          detail="reg" with added plots with regression lines for all batches.
  # NormMoyPool: not used
  # NormInt: not used
  # method: regression method to be used to correct : "linear" or "lowess" or "loess"
  # span: span value for lo(w)ess regression; "none" for linear or default values
  # valNull: to determine what to do with negatively estimated intensities
  # sm_meta: list of information about sample metadata coding
  # min_norm: minimum value accepted for normalisation term (denominator); should be strictly positive
    indfact=which(dimnames(x)[[2]]==fact)
    indtypsamp=which(dimnames(x)[[2]]==sm_meta$sampleType)
    indbatch=which(dimnames(x)[[2]]==sm_meta$batch)
    indinject=which(dimnames(x)[[2]]==sm_meta$injectionOrder)
    lastIon=dim(x)[2]
    indpool=which(x[[sm_meta$sampleType]] %in% sm_meta$sampleTag$pool)# QCpools subscripts in all batches
    valref=apply(as.matrix(x[indpool,(nbid+1):(lastIon)]),2,mean,na.rm=TRUE) # reference value for each ion used to still have the same rought size of values
    nbi=lastIon-nbid # number of ions
    nbb=length(levels(x[[sm_meta$batch]])) # Number of batch(es) = number of levels of factor "batch" (can be =1)
    Xn=data.frame(x[,c(1:nbid)],matrix(0,nrow=nrow(x),ncol=nbi))# initialisation of the corrected dataframe (=initial dataframe)
    dimnames(Xn)=dimnames(x)
    cv=data.frame(matrix(NA,nrow=nbi,ncol=2))# initialisation of dataframe containing CV before and after correction
    dimnames(cv)[[2]]=c("avant","apres")
    if (detail!="reg" && detail!="plot" && detail!="no") {detail="no"}
    pdf(outlog,width=27,height=20)
    cat(nbi," ions ",nbb," batch(es) \n")
    if (detail=="plot") {if(nbb<6){par(mfrow=c(3,3),ask=F,cex=1.5)}else{par(mfrow=c(4,4),ask=F,cex=1.5)}}
    res.ind <- matrix(NA,ncol=nbb,nrow=nbi,dimnames=list(dimnames(x)[[2]][-c(1:nbid)],paste("norm.b",1:nbb,sep="")))
    for (p in 1:nbi) {# for each ion
      labion=dimnames(x)[[2]][p+nbid]
      pools1=x[indpool,p+nbid]
      if(length(which(pools1[!(is.na(pools1))]>0))<2){ # if not enough pools >0 -> no normalisation
        war.note <- paste("Warning: less than 2 pools with values >0 in",labion,"-> no normalisation for this ion.")
        cat(war.note,"\n")
        Xn[,p+nbid] <- x[,p+nbid]
        res.ind[p,] <- rep(0,nbb)
        if (detail=="reg" || detail=="plot" ) {
          par(mfrow=c(2,2),ask=F,cex=1.5)
          plot.new()
          legend("center",war.note)
          minval=min(x[p+nbid],na.rm=TRUE);maxval=max(x[p+nbid],na.rm=TRUE)
          plot( x[[sm_meta$injectionOrder]], x[,p+nbid],col=x[[sm_meta$batch]],ylab=labion,ylim=c(minval,maxval),
          main="No correction",xlab="injection order")
          points(x[[sm_meta$injectionOrder]][indpool],x[indpool,p+nbid],col="maroon",pch=16,cex=1)
        }
      } else {
        if (detail == "reg") {if(nbb<6){par(mfrow=c(3,3),ask=F,cex=1.5)}else{par(mfrow=c(4,4),ask=F,cex=1.5)}}
        if (detail == "plot") {par(mfrow=c(2,2),ask=F,cex=1.5)}
        cv[p,1]=sd(pools1,na.rm=TRUE)/mean(pools1,na.rm=TRUE)# CV before correction
        for (b in 1:nbb) {# for every batch
            indbt = which(x[[sm_meta$batch]]==(levels(x[[sm_meta$batch]])[b])) # subscripts of all samples
            sub=data.frame(x[(x[[sm_meta$batch]]==levels(x[[sm_meta$batch]])[b]),c(indtypsamp,indinject,p+nbid)])
            if (method=="linear") { res.norm = normlinear(sub,detail,valref[p],b,valNull,sm_meta,min_norm)
            } else { if (method=="loess"){ res.norm <- normloess(sub,detail,valref[p],b,span,valNull,sm_meta,min_norm)
              } else { if (method=="lowess"){ res.norm <- normlowess(sub,detail,valref[p],b,span,valNull,sm_meta,min_norm)
                } else {stop("\n--\nNo valid 'method' argument supplied.\nMust be 'linear','loess' or 'lowess'.\n--\n")}
            }}
            Xn[indbt,p+nbid] = res.norm[[1]]
            res.ind[p,b] <- res.norm[[2]]
        }
        # Post correction CV calculation
        pools2=Xn[indpool,p+nbid]
        cv[p,2]=sd(pools2,na.rm=TRUE)/mean(pools2,na.rm=TRUE)
        if (detail=="reg" || detail=="plot" ) {
            # plot before and after correction
            minval=min(cbind(x[p+nbid],Xn[p+nbid]),na.rm=TRUE);maxval=max(cbind(x[p+nbid],Xn[p+nbid]),na.rm=TRUE)
            plot( x[[sm_meta$injectionOrder]], x[,p+nbid],col=x[[sm_meta$batch]],ylab=labion,ylim=c(minval,maxval),
              main=paste0("before correction (CV for pools = ",round(cv[p,1],2),")"),xlab="injection order")
              points(x[[sm_meta$injectionOrder]][indpool],x[indpool,p+nbid],col="maroon",pch=16,cex=1)
            plot(Xn[[sm_meta$injectionOrder]],Xn[,p+nbid],col=x[[sm_meta$batch]],ylab="",ylim=c(minval,maxval),
             main=paste0("after correction (CV for pools = ",round(cv[p,2],2),")"),xlab="injection order")
              points(Xn[[sm_meta$injectionOrder]][indpool],Xn[indpool,p+nbid],col="maroon",pch=16,cex=1)
            suppressWarnings(plot.design( x[c(indtypsamp,indbatch,indfact,p+nbid)],main="factors effect before correction"))
            suppressWarnings(plot.design(Xn[c(indtypsamp,indbatch,indfact,p+nbid)],main="factors effect after correction"))
        }
      }
    }

    if (detail=="reg" || detail=="plot" || detail=="no") {
        if (nbi > 3) {
            # Sum of ions before/after plot
            par(mfrow=c(1,2),ask=F,cex=1.2)
            xsum <- rowSums(x[,(nbid+1):lastIon],na.rm=TRUE)
            Xnsum <- rowSums(Xn[,(nbid+1):lastIon],na.rm=TRUE)
            plot(x[[sm_meta$injectionOrder]],xsum,col=x[[sm_meta$batch]],ylab="sum of variables' intensities",xlab="injection order",
                 ylim=c(min(c(xsum,Xnsum),na.rm=TRUE),max(c(xsum,Xnsum),na.rm=TRUE)),main="Sum of intensities\nBefore correction")
            points(x[[sm_meta$injectionOrder]][indpool],xsum[indpool],col="maroon",pch=16,cex=1.2)
            plot(x[[sm_meta$injectionOrder]],Xnsum,col=x[[sm_meta$batch]],ylab="sum of variables' intensities",xlab="injection order",
                 ylim=c(min(c(xsum,Xnsum),na.rm=TRUE),max(c(xsum,Xnsum),na.rm=TRUE)),main="Sum of intensities\nAfter correction")
            points(x[[sm_meta$injectionOrder]][indpool],Xnsum[indpool],col="maroon",pch=16,cex=1.2)
            # PCA Plot before/after, normed only and ions plot
            par(mfrow=c(3,4),ask=F,cex=1.2) 
            acplight(x[,c(indtypsamp,indbatch,indtypsamp,indfact,(nbid+1):lastIon)],"uv",TRUE)
            norm.ion <- which(colnames(Xn)%in%(rownames(res.ind)[which(rowSums(res.ind)>=1)]))
            acplight(Xn[,c(indtypsamp,indbatch,indtypsamp,indfact,(nbid+1):lastIon)],"uv",TRUE,norm.ion)
            if(length(norm.ion)>0){acplight(Xn[,c(indtypsamp,indbatch,indtypsamp,indfact,norm.ion)],"uv",TRUE)}
            # Before/after boxplot
            par(mfrow=c(1,2),ask=F,cex=1.2) 
            cvplot=cv[!is.na(cv[[1]])&!is.na(cv[[2]]),]
            if(nrow(cvplot)>0){
              boxplot(cvplot[[1]],ylim=c(min(cvplot),max(cvplot)),main="CV of pools before correction")
              boxplot(cvplot[[2]],ylim=c(min(cvplot),max(cvplot)),main="CV of pools after correction")
            }
            dev.off()
        }
    }
  if (nbi<=3) {dev.off()}
  # transposed matrix is return  (format of the initial matrix with ions in rows)
  Xr=Xn[,-c(1:nbid)]; dimnames(Xr)[[1]]=Xn[[1]]
  Xr=t(Xr) ; Xr <- data.frame(ions=rownames(Xr),Xr)
  
  res.norm[[1]] <- Xr ; res.norm[[2]] <- data.frame(metaion,res.ind) ; res.norm[[3]] <- x[,c(1:nbid)]
  names(res.norm) <- c("dataMatrix","variableMetadata","sampleMetadata")
  return(res.norm)
}





acplight <- function(ids, scaling="uv", indiv=FALSE,indcol=NULL) {
    suppressPackageStartupMessages(library(ade4))
    suppressPackageStartupMessages(library(pcaMethods))
  # Make a PCA and plot scores and loadings.
  # First column must contain samples' identifiers.
  # Columns 2 to 4 contain factors to colour the plots. 
  for (i in 1:3) {
    idss <- data.frame(ids)
    idss[,i+1] <- as.character(idss[,i+1])
    idss[which(is.na(idss[,i+1])),i+1] <- "no_modality"
    idss[which(idss[,i+1]=="NA"),i+1] <- "no_modality"
    idss[which(idss[,i+1]==""),i+1] <- "no_modality"
    classe=as.factor(idss[[i+1]])
    idsample=as.character(idss[[1]])
    colour=1:length(levels(classe))
    ions=as.matrix(idss[,5:dim(idss)[2]])
    # Removing ions containing NA (not compatible with standard PCA)
    ions=t(na.omit(t(ions)))
    if(i==1){if(ncol(ions)!=(ncol(idss)-4)){cat("Note:",(ncol(idss)-4)-ncol(ions),"ions were ignored for PCA display due to NA in intensities.\n")}}
    # Scaling choice: "uv","none","pareto"
    object=suppressWarnings(prep(ions, scale=scaling, center=TRUE))
    if(i==1){if(length(which(apply(ions,2,var)==0))>0){cat("Warning: there are",length(which(apply(ions,2,var)==0)),"constant ions.\n")}}
    # ALGO: nipals,svdImpute, Bayesian, svd, probalistic=F
    result <- pca(object, center=F, method="svd", nPcs=2)
    # ADE4 : to plot samples' ellipsoid for each class
    s.class(result@scores, classe, cpoint = 1,xax=1,yax=2,col=colour,sub=sprintf("Scores - PCs %sx%s",1,2), possub="bottomright")
    #s.label(result@loadings,label = ions, cpoint = 0, clabel=0.4, xax=1,yax=2,sub="Loadings",possub="bottomright")
    if(i==1){resulti <- result}
  }
  if(indiv) {
    colour <- rep("darkblue",length(resulti@loadings)) ; if(!is.null(indcol)) {colour[-c(indcol)] <- "red"}
    plot(resulti@loadings,col=colour,main="Loadings",xaxt="n",yaxt="n",pch=20,
         xlab=bquote(PC1-R^2==.(resulti@R2[1])),ylab=bquote(PC2 - R^2 == .(resulti@R2[2])))
    abline(h=0,v=0)}
}



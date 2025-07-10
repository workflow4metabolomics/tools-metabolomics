###############################################################################################################
# batch_correction_3Llauncher                                                                                 #
#                                                                                                             #
# Authors: Jean-Francois MARTIN / Melanie Petera                                                              #
# Starting date: 04-08-2020                                                                                   #
# Based on batch_correction_wrapper.R version 2.91                                                            #
# Version 1: 02-10-2020                                                                                       #
#            - split of tool-linked code and script-linked one                                                #
#            - handling of sample tags' parameters                                                            #
#            - accepting samples beyond pools and samples                                                     #
#            - dealing with special characters in IDs and column names                                        #
#            - adding a min.norm argument to the function                                                     #
#                                                                                                             #
# Input files: dataMatrix.txt, sampleMetadata.txt, variableMetadata.txt (BC only)                             #
# Output files: graph.pdf, corrected table (BC only), diagnostic table (DBC only), variableMetadata (BC only) #
#                                                                                                             #
###############################################################################################################

meth3L <- function(idsample,iddata,sample_type_col_name,injection_order_col_name,batch_col_name,sample_type_tags,
                   factbio,analyse,metaion,detail,method,outlog,span,valnull,
                   rdata_output,dataMatrix_out,variableMetadata_out,out_graph_pdf,out_preNormSummary,min.norm){

## Import function
tab.import <- function(tested.file,tabtype){
  tab.res <- tryCatch(read.table(tested.file,header=TRUE,sep='\t',check.names=FALSE,comment.char = ''), error=conditionMessage)
  if(length(tab.res)==1){
    stop(paste("Could not import the",tabtype,"file. There may be issues in your table integrity.\nCorresponding R error message:\n",tab.res))
  }else{
    tab.comp <- tryCatch(read.table(tested.file,header=TRUE,sep='\t',check.names=FALSE,comment.char = '',quote=""), error=conditionMessage)
    if((length(tab.comp)!=1)&&(dim(tab.res)!=dim(tab.comp))){ # wrong original import due to quotes inside a column name
      return(tab.comp)
    }else{ return(tab.res) }
  }
}

## Reading of input files
idsample=tab.import(idsample,"sampleMetadata")
iddata=tab.import(iddata,"dataMatrix")

### Table match check 
table.check <- match2(iddata,idsample,"sample")
if(length(table.check)>1){check.err(table.check)}

### StockID
samp.id <- stockID(iddata,idsample,"sample")
iddata<-samp.id$dataMatrix ; idsample<-samp.id$Metadata ; samp.id<-samp.id$id.match

### Checking mandatory variables
mand.check <- ""
for(mandcol in c(sample_type_col_name, injection_order_col_name, batch_col_name)){
  if(!(mandcol%in%colnames(idsample))){
    mand.check <- c(mand.check,"\nError: no '",mandcol,"' column in sample metadata.\n",
                    "Note: column names are case-sensitive.\n")
  }
}
if(length(mand.check)>1){
  mand.check <- c(mand.check,"\nFor more information, see the help section or:",
                  "\n http://workflow4metabolomics.org/sites/",
                  "workflow4metabolomics.org/files/files/w4e-2016-data_processing.pdf\n")
  check.err(mand.check)
}

if(analyse == "batch_correction") {
    ## Reading of Metadata Ions file
    metaion=read.table(metaion,header=T,sep='\t',check.names=FALSE,comment.char = '')
    ## Table match check 
    table.check <- c(table.check,match2(iddata,metaion,"variable"))
    ## StockID
    var.id <- stockID(iddata,metaion,"variable")
    iddata<-var.id$dataMatrix ; metaion<-var.id$Metadata ; var.id<-var.id$id.match
}

### Formating
idsample[[1]]=make.names(idsample[[1]])
dimnames(iddata)[[1]]=iddata[[1]]

### Transposition of ions data
idTdata=t(iddata[,2:dim(iddata)[2]])
idTdata=data.frame(dimnames(idTdata)[[1]],idTdata)

### Merge of 2 files (ok even if the two dataframe are not sorted on the same key)
ids=merge(idsample, idTdata, by.x=1, by.y=1)

ids[[batch_col_name]]=as.factor(ids[[batch_col_name]])
nbid=dim(idsample)[2]

### Checking the number of sample and pool

# least 2 samples
if(length(which(ids[[sample_type_col_name]] %in% sample_type_tags$sample))<2){
    table.check <- c(table.check,"\nError: less than 2 samples specified in sample metadata.",
           "\nMake sure this is not due to errors in your ",sample_type_col_name," coding.\n")
}

# least 2 pools per batch for all batchs
B <- rep(0,length(levels(ids[[batch_col_name]])))
for(nbB in 1:length(levels(ids[[batch_col_name]]))){
    B[nbB]<-length(which(ids[which(ids[[batch_col_name]]==(levels(ids[[batch_col_name]])[nbB])),,drop=FALSE][[sample_type_col_name]] %in% sample_type_tags$pool))
}
if(length(which(B>1))==0){
    table.check <- c(table.check,"\nError: less than 2 pools specified in at least one batch in sample metadata.",
           "\nMake sure this is not due to errors in your ",sample_type_col_name," coding.\n")
}

### Checking the unicity of samples and variables
uni.check <- function(tested.tab,tabtype,err.obj){
  unicity <- duplicated(tested.tab[,1])
  if(sum(unicity)>0){
    #Sending back an explicit error
    duptable <- t(t(table(tested.tab[,1][unicity])+1))
    err.obj <- c(err.obj,paste0("\n-------\nError: your '",tabtype,"' IDs contain duplicates:\n"),
                 paste(rownames(duptable),duptable,sep=": ",collapse="\n"),
                 "\nSince identifiers are meant to be unique, please check your data.\n-------\n")
  }
  return(err.obj)
}
table.check <- uni.check(ids,"sample",table.check)
if(analyse == "batch_correction"){table.check <- uni.check(metaion,"variable",table.check)}

## error check
check.err(table.check)


### BC/DBC-specific processing

# Gathering mandatory information in a single object
sm.meta <- list(batch=batch_col_name, injectionOrder=injection_order_col_name, sampleType=sample_type_col_name, sampleTag=sample_type_tags)

if(analyse == "batch_correction") {
    ## Launch
    res = norm_QCpool(ids,nbid,outlog,factbio,metaion,detail,FALSE,FALSE,method,span,valnull,sm.meta,min.norm)
    ## Get back original IDs
    var.id <- reproduceID(res[[1]],res[[2]],"variable",var.id)
    res[[1]] <- var.id$dataMatrix ; res[[2]] <- var.id$Metadata
    samp.id <- reproduceID(res[[1]],res[[3]],"sample",samp.id)
    res[[1]] <- samp.id$dataMatrix ; res[[3]] <- samp.id$Metadata
    ## Save files
    save(res, file=rdata_output)
    write.table(res[[1]], file=dataMatrix_out, sep = '\t', row.names=FALSE, quote=FALSE)
    write.table(res[[2]], file=variableMetadata_out, sep = '\t', row.names=FALSE, quote=FALSE)
}else{
    ## Launch
    plotsituation(ids,nbid,out_graph_pdf,out_preNormSummary,factbio,span,sm.meta)
}

}#end of meth3L
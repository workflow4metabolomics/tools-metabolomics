################################################################################################
# TABLE MERGE                                                                                  #
#                                                                                              #
# User: Galaxy                                                                                 #
# Starting date: 16-04-2015                                                                    #
# V-0.1: First version of merge code                                                           #
# V-0.2: Addition of data check and handling of special characters                             #
#                                                                                              #
#                                                                                              #
# Input files: dataMatrix ; Metadata file                                                      #
# Output files: dataMatrix ; Metadata file                                                     #
#                                                                                              #
# Dependencies: RcheckLibrary.R ; miniTools.R                                                  #
#                                                                                              #
################################################################################################

# Parameters (for dev)
if(FALSE){
  DM.name <- "dataMatrix_CleanIons_CleanEch.txt"
  meta.name <- "sampleMetadata_CleanEch.txt"
  metype <- "sample"
  output <- "Combined_${Metadata_in.name}"
}



tab.merge <- function(DM.name,meta.name,metype,output){
  # This function allows to merge the dataMatrix with one metadata table.
  #
  # Parameters:
  # - DM.name, meta.name: dataMatrix and metadata files' access respectively
  # - metype: "sample" or "variable" depending on metadata content
  # - output: output file's access
  
  
# Input --------------------------------------------------------------

DM <- read.table(DM.name,header=TRUE,sep="\t",check.names=FALSE)
meta <- read.table(meta.name,header=TRUE,sep="\t",check.names=FALSE,colClasses="character")

# Table match check 
table.check <- match2(DM,meta,metype)
check.err(table.check)

# StockID
meta.id <- stockID(DM,meta,metype)
DM<-meta.id$dataMatrix ; meta<-meta.id$Metadata ; meta.id<-meta.id$id.match


# Merging tables -----------------------------------------------------

if(metype=="sample"){
  ori.DM <- DM
  rownames(DM) <- DM[,1]
  DM <- DM[,-1]
  DM <- t(DM)
  DM <- data.frame(sample=row.names(DM),DM,check.names=FALSE)
  rownames(DM) <- NULL
}

comb.data <- merge(x=meta,y=DM,by.x=1,by.y=1)


# Output -------------------------------------------------------------

# Getting back original identifiers
if(metype=="sample"){
  id.ori <- reproduceID(ori.DM,comb.data,metype,meta.id)
}else{
  id.ori <- reproduceID(DM,comb.data,metype,meta.id)
}
comb.data <- id.ori$Metadata

# Writing the table
write.table(comb.data,output,sep="\t",quote=FALSE,row.names=FALSE)


} # End of tab.merge


# Typical function call
# tab.merge(DM.name,meta.name,metype,output)

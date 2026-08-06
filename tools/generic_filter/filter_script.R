################################################################################################
# GENERIC FILTERS                                                                              #
#                                                                                              #
# User: Galaxy                                                                                 #
# Starting date: 03-09-2014                                                                    #
# V-1.0: Restriction of old filter script to Filter according to factors                       #
# V-1.1: Choice of metadata table for filtering added ; data check added ; handling of NA ;    #
#        check for minimum remaining data                                                      #
# V-1.2: Minor modifications in script layout                                                  #
# V-2.0: Addition of numerical filter                                                          #
# V-2.1: Handling special characters                                                           #
#                                                                                              #
#                                                                                              #
# Input files: dataMatrix ; sampleMetadata ; variableMetadata                                  #
# Output files: dataMatrix ; sampleMetadata ; variableMetadata                                 #
#                                                                                              #
################################################################################################

# Parameters (for dev)
if(FALSE){
  
  ion.file.in <- "test/ressources/inputs/ex_data_IONS.txt"  #tab file
  meta.samp.file.in <- "test/ressources/inputs/ex_data_PROTOCOLE1.txt"  #tab file
  meta.ion.file.in <- "test/ressources/inputs/ex_data_METAION.txt"  #tab file
  
  ion.file.out <- "test/ressources/outputs/ex_data_IONS_fl.txt"  #tab file
  meta.samp.file.out <- "test/ressources/outputs/ex_data_PROTOCOLE1_fl.txt"  #tab file
  meta.ion.file.out <- "test/ressources/outputs/ex_data_METAION_fl.txt"  #tab file

NUM <- TRUE ; if(NUM){ls.num<-list(c("sample","injectionOrder","upper","20"),c("variable","var1","extremity","0.12","500"))}else{ls.num<-NULL}
  
FACT <- TRUE ; if(FACT){ls.fact<-list(c("centre","C","sample"),c("var2","A","variable"))}else{ls.fact<-NULL}
  
}

filters <- function(ion.file.in, meta.samp.file.in, meta.ion.file.in,
                    NUM, ls.num, FACT, ls.fact,
                    ion.file.out, meta.samp.file.out, meta.ion.file.out){
  # This function allows to filter variables and samples according to factors or numerical values. 
  # It needs 3 datasets: the data matrix, the variables' metadata, the samples' metadata. 
  # It generates 3 new datasets corresponding to the 3 inputs filtered. 
  #
  # Parameters:
  # - xxx.in: input files' access
  # - xxx.out: output files' access
  # - NUM: filter according to numerical variables yes/no
  # | > ls.num: numerical variables' list for filter
  # - FACT: filter according to factors yes/no
  # | > ls.fact: factors' list for filter
  
  
# Input -----------------------------------------------------------------------------------

ion.data <- read.table(ion.file.in, sep = "\t", header = TRUE, check.names = FALSE)
meta.samp.data <- read.table(meta.samp.file.in, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = TRUE)
meta.ion.data <- read.table(meta.ion.file.in, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = TRUE)

# Error vector
err.stock <- "\n"


# Table match check 
table.check <- match3(ion.data,meta.samp.data,meta.ion.data)
check_err(table.check)

# StockID
samp.id <- stock_id(ion.data,meta.samp.data,"sample")
ion.data <- samp.id$dataMatrix
meta.samp.data <- samp.id$Metadata
samp.id <- samp.id$id.match



# Function 1: Filter according to numerical variables -------------------------------------
# Allows to delete all elements corresponding to defined values of designated variables.
if(NUM){
  
  # For each numerical variable to filter
  for(i in 1:length(ls.num)){
    
    # Which metadata table is concerned
    if(ls.num[[i]][1]=="sample"){metadata <- meta.samp.data}else{metadata <- meta.ion.data}
    
    # Checking the columns and factors variables
    numcol <- which(colnames(metadata)==ls.num[[i]][2])
    if(length(numcol)==0) {
      err.stock <- c(err.stock,"\n-------",
                     "\nWarning: no '",ls.num[[i]][2],"' column detected in ",ls.num[[i]][1],
                     " metadata!","\nFiltering impossible for this variable.\n-------\n") 
    }else{
      if(!is.numeric(metadata[,numcol])){
        err.stock <- c(err.stock,"\n-------",
                       "\nWarning: column '",ls.num[[i]][2],"' in ",ls.num[[i]][1],
                       " metadata is not a numerical variable!",
                       "\nNumerical filtering impossible for this variable.\n-------\n")
      }else{
        
        # Filtering
        if(ls.num[[i]][3]=="lower"){
          toremove <- which(metadata[,numcol]<as.numeric(ls.num[[i]][4]))
          if(length(toremove)!=0){
            metadata <- metadata[-c(toremove),]
          }
        }else{if(ls.num[[i]][3]=="upper"){
          toremove <- which(metadata[,numcol]>as.numeric(ls.num[[i]][4]))
          if(length(toremove)!=0){
            metadata <- metadata[-c(toremove),]
          }
        }else{if(ls.num[[i]][3]=="between"){
          toremove <- (metadata[,numcol]>as.numeric(ls.num[[i]][4]))+(metadata[,numcol]<as.numeric(ls.num[[i]][5]))
          toremove <- which(toremove==2)
          if(length(toremove)!=0){
            metadata <- metadata[-c(toremove),]
          }
        }else{if(ls.num[[i]][3]=="extremity"){
          toremove <- c(which(metadata[,numcol]<as.numeric(ls.num[[i]][4])),
                        which(metadata[,numcol]>as.numeric(ls.num[[i]][5])))
          if(length(toremove)!=0){
            metadata <- metadata[-c(toremove),]
          }
        }}}}
        
        # Extension to the tables
        if(ls.num[[i]][1]=="sample"){
          meta.samp.data <- metadata
          ion.data <- ion.data[,c(1,which(colnames(ion.data)%in%meta.samp.data[,1]))]
        }else{
          meta.ion.data <- metadata
          ion.data <- ion.data[which(ion.data[,1]%in%meta.ion.data[,1]),]
        }
        
      }}}
  
} # end if(NUM)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - -



# Function 2: Filter according to factors -------------------------------------------------
# Allows to delete all elements corresponding to selected value of designated factor.
if(FACT){

  # For each factor to filter
  for(i in 1:length(ls.fact)){
    
	# Which metadata table is concerned
	if(ls.fact[[i]][3]=="sample"){metadata <- meta.samp.data}else{metadata <- meta.ion.data}
	
    # Checking the columns and factors variables
    numcol <- which(colnames(metadata)==ls.fact[[i]][1])
    if(length(numcol)==0) {
    err.stock <- c(err.stock,"\n-------",
                   "\nWarning: no '",ls.fact[[i]][1],"' column detected in ",ls.fact[[i]][3],
                   " metadata!","\nFiltering impossible for this factor.\n-------\n") 
    }else{
    if((!(ls.fact[[i]][2]%in%levels(as.factor(metadata[,numcol]))))&((ls.fact[[i]][2]!="NA")|(length(which(is.na(metadata[,numcol])))==0))){
      err.stock <- c(err.stock,"\n-------",
                     "\nWarning: no '",ls.fact[[i]][2],"' level detected in '",
                     ls.fact[[i]][1],"' column (",ls.fact[[i]][3]," metadata)!\n",
					 "Filtering impossible for this factor.\n-------\n")
    }else{
      
    # Filtering
    if(length(which(metadata[,numcol]==ls.fact[[i]][2]))!=0){ #if the level still exists in the data
      metadata <- metadata[-c(which(metadata[,numcol]==ls.fact[[i]][2])),]
	}else{ #to treat the special case of "NA" level
	  if(ls.fact[[i]][2]=="NA"){metadata <- metadata[-c(which(is.na(metadata[,numcol]))),]}
	}
	
	# Extension to the tables
	if(ls.fact[[i]][3]=="sample"){
	  meta.samp.data <- metadata
      ion.data <- ion.data[,c(1,which(colnames(ion.data)%in%meta.samp.data[,1]))]
	}else{
	  meta.ion.data <- metadata
      ion.data <- ion.data[which(ion.data[,1]%in%meta.ion.data[,1]),]
	}

  }}}

} # end if(FACT)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - -




# Check if at least one sample and one variable remain ------------------------------------

if(nrow(meta.samp.data)==0){
  stop("\n /!\\ Your filtering options lead to no more sample in your data matrix!\n",
       "Think about reducing your number of filter.")
}

if(nrow(meta.ion.data)==0){
  stop("\n /!\\ Your filtering options lead to no more variable in your data matrix!\n",
       "Think about reducing your number of filter.")
}

# Output ----------------------------------------------------------------------------------

# Getting back original identifiers
id.ori <- reproduce_id(ion.data,meta.samp.data,"sample",samp.id)
ion.data <- id.ori$dataMatrix
meta.samp.data <- id.ori$Metadata


# Error checking
if(length(err.stock)>1){
  stop(err.stock)
}else{

write.table(ion.data, ion.file.out, sep="\t", row.names=FALSE, quote=FALSE)
write.table(meta.samp.data, meta.samp.file.out, sep="\t", row.names=FALSE, quote=FALSE)
write.table(meta.ion.data, meta.ion.file.out, sep="\t", row.names=FALSE, quote=FALSE)

}


} # end of filters function


# Typical function call
#filters(ion.file.in, meta.samp.file.in, meta.ion.file.in, 
#        NUM, ls.num, FACT, ls.fact,
#        ion.file.out, meta.samp.file.out, meta.ion.file.out)


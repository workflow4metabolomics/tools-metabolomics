######################################################
# R check library
# Coded by: M.Petera, 
# - -
# R functions to use in R scripts  
# (management of various generic subroutines)
# - -
# V0: script structure + first functions
# V1: More detailed error messages in match functions
######################################################


# Generic function to return an error if problems have been encountered - - - -

check.err <- function(err.stock){
	
	# err.stock = vector of results returned by check functions
	
	if(length(err.stock)!=0){ stop("\n- - - - - - - - -\n",err.stock,"\n- - - - - - - - -\n") }
	
}




# Table match check functions - - - - - - - - - - - - - - - - - - - - - - - - -

# To check if dataMatrix and (variable or sample)Metadata match regarding identifiers
match2 <- function(dataMatrix, Metadata, Mtype){

	# dataMatrix = data.frame containing dataMatrix
	# Metadata = data.frame containing sampleMetadata or variableMetadata
	# Mtype = "sample" or "variable" depending on Metadata content
	
	err.stock <- NULL # error vector

	id2 <- Metadata[,1]
	if(Mtype=="sample"){ id1 <- colnames(dataMatrix)[-1] }
	if(Mtype=="variable"){ id1 <- dataMatrix[,1] }
	
	if( length(which(id1%in%id2))!=length(id1) || length(which(id2%in%id1))!=length(id2) ){
		err.stock <- c("\nData matrix and ",Mtype," metadata do not match regarding ",Mtype," identifiers.")
		if(length(which(id1%in%id2))!=length(id1)){
			if(length(which(!(id1%in%id2)))<4){ err.stock <- c(err.stock,"\n    The ")
			}else{ err.stock <- c(err.stock,"\n    For example, the ") }
			err.stock <- c(err.stock,"following identifiers found in the data matrix\n",
							"    do not appear in the ",Mtype," metadata file:\n")
			identif <- id1[which(!(id1%in%id2))][1:min(3,length(which(!(id1%in%id2))))]
			err.stock <- c(err.stock,"    ",paste(identif,collapse="\n    "),"\n")
		}
		if(length(which(id2%in%id1))!=length(id2)){
			if(length(which(!(id2%in%id1)))<4){ err.stock <- c(err.stock,"\n    The ")
			}else{ err.stock <- c(err.stock,"\n    For example, the ") }
			err.stock <- c(err.stock,"following identifiers found in the ",Mtype," metadata file\n",
							"    do not appear in the data matrix:\n")
			identif <- id2[which(!(id2%in%id1))][1:min(3,length(which(!(id2%in%id1))))]
			err.stock <- c(err.stock,"    ",paste(identif,collapse="\n    "),"\n")
		}
		err.stock <- c(err.stock,"\nPlease check your data.\n")
	}
	
	return(err.stock)
	
}

# To check if the 3 standard tables match regarding identifiers
match3 <- function(dataMatrix, sampleMetadata, variableMetadata){
	
	# dataMatrix = data.frame containing dataMatrix
	# sampleMetadata = data.frame containing sampleMetadata
	# variableMetadata = data.frame containing variableMetadata
	
	err.stock <- NULL # error vector
	
	id1 <- colnames(dataMatrix)[-1]
	id2 <- sampleMetadata[,1]
	id3 <- dataMatrix[,1]
	id4 <- variableMetadata[,1]
	
	if( length(which(id1%in%id2))!=length(id1) || length(which(id2%in%id1))!=length(id2) ){
		err.stock <- c(err.stock,"\nData matrix and sample metadata do not match regarding sample identifiers.")
		if(length(which(id1%in%id2))!=length(id1)){
			if(length(which(!(id1%in%id2)))<4){ err.stock <- c(err.stock,"\n    The ")
			}else{ err.stock <- c(err.stock,"\n    For example, the ") }
			err.stock <- c(err.stock,"following identifiers found in the data matrix\n",
							"    do not appear in the sample metadata file:\n")
			identif <- id1[which(!(id1%in%id2))][1:min(3,length(which(!(id1%in%id2))))]
			err.stock <- c(err.stock,"    ",paste(identif,collapse="\n    "),"\n")
		}
		if(length(which(id2%in%id1))!=length(id2)){
			if(length(which(!(id2%in%id1)))<4){ err.stock <- c(err.stock,"\n    The ")
			}else{ err.stock <- c(err.stock,"\n    For example, the ") }
			err.stock <- c(err.stock,"following identifiers found in the sample metadata file\n",
							"    do not appear in the data matrix:\n")
			identif <- id2[which(!(id2%in%id1))][1:min(3,length(which(!(id2%in%id1))))]
			err.stock <- c(err.stock,"    ",paste(identif,collapse="\n    "),"\n")
		}
	}
	
	if( length(which(id3%in%id4))!=length(id3) || length(which(id4%in%id3))!=length(id4) ){
		err.stock <- c(err.stock,"\nData matrix and variable metadata do not match regarding variable identifiers.")
		if(length(which(id3%in%id4))!=length(id3)){
			if(length(which(!(id3%in%id4)))<4){ err.stock <- c(err.stock,"\n    The ")
			}else{ err.stock <- c(err.stock,"\n    For example, the ") }
			err.stock <- c(err.stock,"following identifiers found in the data matrix\n",
							"    do not appear in the variable metadata file:\n")
			identif <- id3[which(!(id3%in%id4))][1:min(3,length(which(!(id3%in%id4))))]
			err.stock <- c(err.stock,"    ",paste(identif,collapse="\n    "),"\n")
		}
		if(length(which(id4%in%id3))!=length(id4)){
			if(length(which(!(id4%in%id3)))<4){ err.stock <- c(err.stock,"\n    The ")
			}else{ err.stock <- c(err.stock,"\n    For example, the ") }
			err.stock <- c(err.stock,"following identifiers found in the variable metadata file\n",
							"    do not appear in the data matrix:\n")
			identif <- id4[which(!(id4%in%id3))][1:min(3,length(which(!(id4%in%id3))))]
			err.stock <- c(err.stock,"    ",paste(identif,collapse="\n    "),"\n")
		}
	}
	
	if(length(err.stock)!=0){ err.stock <- c(err.stock,"\nPlease check your data.\n") }
	
	return(err.stock)
	
}
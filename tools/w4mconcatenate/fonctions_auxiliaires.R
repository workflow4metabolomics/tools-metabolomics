#--------------------------------------------------------------------------------------------------------------------------------------------------------------
check_features <- function(M1, M2) {
  #M1/M2 = data.frame containing sampleMetadata or variableMetadata
  #check the variables in the 2 matrices .
  #returns the names of the columns in the two metadata
  
  colnames_1 <- colnames(M1)
  colnames_2 <- colnames(M2)
  samecolumns <- intersect(colnames_1, colnames_2)
  
  if (is.null(samecolumns)) {
    cat("\nWarning: There are no features in common \n")
  }
  
  return(samecolumns)
}

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

compare_same_columns <- function(M1, M2, same, choice_keep, Keep, tab1, tab2) {
  
  #M1/M2 = data.frame containing sampleMetadata or variableMetadata
  #same = list of column names with the same name in M1 and M2
  #choice_keep = choice of keeping columns with the same no or keeping just one
  #keep = keep the column in M1 or M2
  #tab1/tab2 = Suffix for Metadata 1/2
  
  #Check that the variables in the 2 matrices have the same values and 
  #If not, they are renamed or deleted according to the parameters chosen by the user.
  #returns the two modified metadata 
  
  compare_results <- list()
  non_identical_columns_v <- c()
  
  
  #Creation of 2 sub-tables with shared individuals and variables
  common_individuals <- intersect(M1[, 1] , M2[, 1])
  common_columns <- intersect(colnames(M1), colnames(M2))
  
  M1_common <- subset(M1, M1[, 1] %in% common_individuals)
  M2_common <- subset(M2, M2[, 1] %in% common_individuals)
  
  
  M1_common <- M1_common[, common_columns]
  M2_common <- M2_common[, common_columns]
  
  

  common_columns <- common_columns[-1]# delete id column
  
  
  for (col_name in common_columns) {
    #Check that the columns are identical, then delete them from M2
    if (!identical(M1_common[[col_name]], M2_common[[col_name]])) {
      
      non_identical_columns_v <- c(non_identical_columns_v, col_name)
    
      #otherwise store the columns where the values are not the same in non_identical_columns
    } else { 
      M2 <- M2[, -which(colnames(M2) == col_name)]
    }
    
  }

  #if the list of columns that do not take the same values is null, we return M1/M2
  if (is.null(non_identical_columns_v)) {
    
    result <- list(M1 = M1, M2 = M2)
    
    return(result)
    
    
  } else {
     
      for (non_identical_columns in non_identical_columns_v) {
        
      #If we decide to keep the 2 columns and they do not take the same values, we change their names by adding a suffix.
        
      if (choice_keep == "yes") {#keep both columns and give them a new name
         
        new_name <- paste(tab1, non_identical_columns, sep = "_")
        colnames(M1)[colnames(M1) == non_identical_columns] <- new_name
        
        new_name <- paste(tab2, non_identical_columns, sep = "_")
        colnames(M2)[colnames(M2) == non_identical_columns] <- new_name
        
      }
      
      if (choice_keep == "no") {#Keep only one and delete the other
        if (Keep == 1) {
          M2 <- M2[, -which(colnames(M2) == non_identical_columns)]
          new_name <- paste(tab1, non_identical_columns, sep = "_")
          colnames(M1)[colnames(M1) == non_identical_columns] <- new_name
          
        }
        if (Keep == 2) {
          M1 <- M1[, -which(colnames(M1) == non_identical_columns)]
          new_name <- paste(tab2, non_identical_columns, sep = "_")
          colnames(M2)[colnames(M2) == non_identical_columns] <- new_name
        }
        
      }
      } 
    }
  

  result <- list(M1 = M1, M2 = M2)
  return(result)
  
}
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
match2_bis <- function(Metadata_1, Metadata_2, Mtype) {
  
  
  #Metadata1/Metadata2 = data.frame containing sampleMetadata or variableMetadata
  #Mtype = "sample" or "variable" depending on Metadata content
  #To check if metadata1 and metadata2  match regarding identifiers
  #returns a vector containing an error message if the identifiers are not all the same in the two metadatas
  err.stock <- NULL#error vector
  
  
  id2 <- Metadata_1[, 1]
  id1 <- Metadata_2[, 1] 
  
  if (length(which(id1 %in% id2)) != length(id1) || length(which(id2 %in% id1)) != length(id2)) {
    err.stock <- c("\n", Mtype, "Metadata_1 and ", Mtype, "Metadata_2 do not match regarding Metadata_2 identifiers.")
    if (length(which(id1 %in% id2)) != length(id1)) {
      if (length(which( ! (id1 %in% id2))) < 4) {
        err.stock <- c(err.stock, "\n    The ")
      } else {
        err.stock <- c(err.stock, "\n    For example, the ") }
      
        err.stock <- c(err.stock, "following identifiers found in the ", Mtype, "Metadata_1 file\n",
                     "    do not appear in the ", Mtype, " Metadata_2 file:\n")
        identif <- id1[which( ! (id1 %in% id2))][1 : min(3, length(which( ! (id1 %in% id2))))]
        err.stock <- c(err.stock, "    ", paste(identif, collapse = "\n    "), "\n")
    }
    if (length(which(id2 %in% id1)) != length(id2)){
      if (length(which( ! (id2 %in% id1))) < 4) { 
        err.stock <- c(err.stock, "\n    The ")
      } else { 
        err.stock <- c(err.stock, "\n    For example, the ") }
      err.stock <- c(err.stock, "following identifiers found in the ", Mtype, " Metadata_2 file\n",
                     "    do not appear in the", Mtype, " Metadata_1 file:\n")
      identif <- id2[which( ! (id2 %in% id1))][1 : min(3, length(which( ! (id2 %in% id1))))]
      err.stock <- c(err.stock, "    ", paste(identif, collapse = "\n    "), "\n")
    }
    err.stock <- c(err.stock, "\nPlease check your data.\n")
  }
  
  return(err.stock)
  
}
rm(list=ls())
#setwd("Y:/Données")


# Chargement des library et des ressources 

library ("W4MRUtils")
library(dplyr)


source_local("concatenation.R")
source_local("fonctions_auxiliaires.R")


para <- W4MRUtils::parse_args(args = commandArgs()) 

#para <- list(dataMatrix_1 = "Input_Unique_Test_1-2-3-4-5_DM1.txt", dataMatrix_2 = "Input_Unique_Test_1-2-3_DM2.txt", metadata_1 = "Input_Unique_Test_1_M1.txt", metadata_2 = "Input_Unique_Test_1_M2.txt", 
#            type = "sample", tab1 = "Tab1", tab2 = "Tab2", concatenation = "unique", choice_keep = "oui", keep = 0, dataMatrix_1_out = "Datamatrix_1.file.tsv", dataMatrix_2_out = "Datamatrix_2.file.tsv",
#            metadata_out = "Metadata_concatenate.tsv")


cat('\nJob starting time:\n',format(Sys.time(), "%a %d %b %Y %X"),
    '\n\n--------------------------------------------------------------------', 
    '\nParameters used in "concatenation":\n\n')
print(para)
cat('--------------------------------------------------------------------\n\n')




#Lancement de l'outil

A <- W4MRUtils::import2(para$dataMatrix_1, para$metadata_1, para$type, disable_comm = FALSE)
B <- W4MRUtils::import2(para$dataMatrix_2, para$metadata_2, para$type, disable_comm = FALSE)  

DM1 <- A$dataMatrix
M1 <- A$metadata

DM2 <- B$dataMatrix
M2 <- B$metadata

#para = list(dataMatrix_1 = "Input_Unique_Test_1-2-3-4-5_DM1.txt" , dataMatrix_2 = "Input_Unique_Test_1-2-3_DM2.txt", Metadata_1 = "Input_Unique_Test_1_M1.txt", Metadata_2 = "Input_Unique_Test_1_M1.txt", 
#            type = "sample", tab1 = "Tab1", tab2= "Tab2", concatenations = "unique", choice_keep = "oui", Keep = 0)






result_tables <- concat(DM1, M1, DM2, M2, para$type, para$tab1, 
              para$tab2, para$concatenation, para$choice_keep, para$keep)

write.table(result_tables[[1]], file = para$metadata_out, sep="\t", row.names = FALSE, quote = FALSE)
write.table(result_tables[[2]], file = para$dataMatrix_1_out, sep="\t", row.names = FALSE, quote = FALSE)
write.table(result_tables[[3]], file = para$dataMatrix_2_out, sep="\t", row.names = FALSE, quote = FALSE)

cat('\n--------------------------------------------------------------------',
    '\nInformation about R (version, Operating System, attached or loaded packages):\n\n')
sessionInfo()
cat('--------------------------------------------------------------------\n',
    '\nJob ending time:\n',format(Sys.time(), "%a %d %b %Y %X"))




























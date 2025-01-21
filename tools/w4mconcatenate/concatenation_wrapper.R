rm(list = ls())
# Chargement des library et des ressources

library("W4MRUtils")
library(dplyr)


source_local("concatenation.R")
source_local("fonctions_auxiliaires.R")


para <- W4MRUtils::parse_args(args = commandArgs())


cat(
    "\nJob starting time:\n", format(Sys.time(), "%a %d %b %Y %X"),
    "\n\n--------------------------------------------------------------------",
    "\nParameters used by the 'W4M concatenate' tool:\n\n"
)
print(para)
cat("--------------------------------------------------------------------\n\n")


# Lancement de l'outil

A <- W4MRUtils::import2(para$dataMatrix_1, para$metadata_1, para$type, disable_comm = FALSE)
B <- W4MRUtils::import2(para$dataMatrix_2, para$metadata_2, para$type, disable_comm = FALSE)

DM1 <- A$dataMatrix
M1 <- A$metadata

DM2 <- B$dataMatrix
M2 <- B$metadata


result_tables <- concat(DM1, M1, DM2, M2, para$type, para$tab1, para$tab2, para$concatenation, para$choice_keep, para$keep)

write.table(result_tables[[1]], file = para$metadata_out, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(result_tables[[2]], file = para$dataMatrix_1_out, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(result_tables[[3]], file = para$dataMatrix_2_out, sep = "\t", row.names = FALSE, quote = FALSE)

cat(
    "\n--------------------------------------------------------------------",
    "\nInformation about R (version, Operating System, attached or loaded packages):\n\n"
)
sessionInfo()
cat(
    "--------------------------------------------------------------------\n",
    "\nJob ending time:\n", format(Sys.time(), "%a %d %b %Y %X")
)

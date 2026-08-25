rm(list = ls())
# Chargement des library et des ressources

library(OptiLCMS)

para <- W4MRUtils::parse_args(args = commandArgs())

cat(
  "\nJob starting time:\n", format(Sys.time(), "%a %d %b %Y %X"),
  "\n\n--------------------------------------------------------------------",
  "\nParameters used by the 'W4M concatenate' tool:\n\n"
)
print(para)
cat("--------------------------------------------------------------------\n\n")

# Lancement de l'outil

# Extraction des ROIs depuis les échantillons
DataFiles <- list.files(para$image, full.names = TRUE)
mSet <- PerformROIExtraction(datapath = DataFiles, rt.idx = para$rt_idx, rmConts = para$rmConts)

# Optimisation des paramètres en fonction des ROIs extraits
best_params <- PerformParamsOptimization(mSet, param = SetPeakParam(platform = para$platform), ncore = para$ncore)

# Génération d'une visualisation PNG de mSet
png(filename = "chromatogram.png", width = 800, height = 600)
PlotROIs(mSet)
dev.off()

# Sauvegarde des paramètres optimisés dans un fichier texte
write.table(
  best_params,
  file = "optimized_parameters.txt",
  sep = "\t",
  col.names = NA,
  quote = FALSE
)

cat(
  "\n--------------------------------------------------------------------",
  "\nInformation about R (version, Operating System, attached or loaded packages):\n\n"
)
sessionInfo()
cat(
  "--------------------------------------------------------------------\n",
  "\nJob ending time:\n", format(Sys.time(), "%a %d %b %Y %X")
)

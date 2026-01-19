#' Summary plot for cross validation.
#' 
#' Creates a summary ploy by comparing predictive abilities for each phenotype and each model
#'
#' @param predictionsKFoldsTable Table with information from cross validation: GEBVs for each fold and each repetition.
#' @return None
#' @import ggplot2
#' @import dplyr 
#'
#' @export
#-------------------------------------------------------------
#--- Create output tables and boxplots for best model for each trait
#-------------------------------------------------------------
gs_summaryPlot <- function (predictionsKFoldsTable) {
	gs_plotsHCL (predictionsKFoldsTable, "Predictive_ability", "GS Summary plot")
}

#-------------------------------------------------------------
#
#-------------------------------------------------------------
gs_outputs <- function (resultsDir="", kfoldsTable=""){
	if (resultsDir!="")
		kfoldsTable  = createSummaryTables (resultsDir)
	gs_summaryPlot (kfoldsTable)

}

#-------------------------------------------------------------
#--- Create output tables and boxplots for best model for each trait
#-------------------------------------------------------------
gs_plotsHCL <- function (HCLTable, columnName, TITLE) {
	nBaseTraits = length (unique (HCLTable [,1]))
	message (">>> Creating output plots and tables for best model ...")

	# Summary plot boxplot best models for traits
	meansTable = HCLTable %>% dplyr::group_by (Traits, Models) %>% 
		summarize (Means=mean(get(columnName)), .groups="keep") 

	# Get best model for each trait
	bestTable = meansTable %>% group_by (Traits) %>% slice_max (Means) %>% slice_head

	# Create table with correlations from best models
	getBestCorr <- function (trait, model) 
		return (dplyr::filter (HCLTable, Traits==trait, Models==model))
	bestCorrs = do.call (rbind, mapply (getBestCorr, bestTable$Traits, bestTable$Models, SIMPLIFY=F))

	traitColors = unlist (lapply(1:nBaseTraits, function(x) rep(x,3))) 
	nTraits     = length (unique (bestCorrs$Traits))
	traitColors = traitColors [1: nTraits] 
	outPlotname   = "out-CrossValidation-ModelsTraits.pdf"
	ggplot (bestCorrs, aes(x=HCL_component, y=get(columnName))) + ylim (0,1) +
		    geom_boxplot (alpha=0.3, fill=traitColors) + 
			theme(axis.text.x=element_text(angle=0, hjust=1)) + labs (title="Best GS models for color traits in tetraploid potato") + 
			stat_summary(geom='text', label=round (bestTable$Means,2) , fun=max, vjust = -1, size=3, col="blue") +
			stat_summary(geom='text', label=bestTable$Models, fun=min, vjust = 1.5, angle=0,size=3, col="blue") + 
			facet_wrap (~Prefix, ncol=nBaseTraits) 
	ggsave (outPlotname, width=11)
}
#-------------------------------------------------------------
#-------------------------------------------------------------
# Create summary tables from GS results of each phenotype
# Create a GEBVs table and kfolds table, the last is returned to be graphed
#-------------------------------------------------------------
createSummaryTables <- function (inputDir) {
	message (">>> Creating table from cross validations predictive abilities...")

	fileCrossValPredictions  = "out-CrossValidation-kfolds.csv"
	fileGEBVsPredictions     = "out-GEBVs-BestModel-TARGET.csv"

	predictionsKFoldsTable   = data.frame ()

#	phenoNames               = colnames (read.csv (phenotypeFile)[,-1])
	phenoNames               = list.files (inputDir)

	predictionsTableList     = list()
	for (pheno in phenoNames) {
		traitPrefix = strsplit (pheno, "[.]")[[1]][1]
		traitHCL    = strsplit (pheno, "[.]")[[1]][2]

		# Construct table of cross validation means, checking if CV data exists
		predictionsFilePath = sprintf ("%s/%s/%s", inputDir, pheno, fileCrossValPredictions)
		if (!file.exists (predictionsFilePath)) {
			message ("WARNING: file ", predictionsFilePath, " not found!")
			cvKFoldsTable  = as.data.frame (list (Models="",Predictive_ability=0))
		}else {
			cvKFoldsTable  = read.csv (predictionsFilePath)

			# Build GEBVs table for all phenotypes
			filePathGEBVsPredictions = sprintf ("%s/%s/%s", inputDir, pheno, fileGEBVsPredictions)
			predictionsTable         = read.csv (filePathGEBVsPredictions, row.names=1)
			names (predictionsTable) = c (pheno, sprintf ("%s.GEBV", pheno))
			predictionsTableList     = append (predictionsTableList, predictionsTable)
		}
		
		# Add prefix and component 
		dfModels = data.frame (Prefix=traitPrefix, HCL_component=traitHCL, Traits=pheno, cvKFoldsTable)
		predictionsKFoldsTable = rbind (predictionsKFoldsTable, dfModels)
	}

	# Write GEBVs file
	#outFile = gsub (".csv", "-ALL.csv", fileGEBVsPredictions)
	outFile = addLabel (fileGEBVsPredictions, "ALL")
	df      = data.frame (SAMPLES=rownames(predictionsTable), predictionsTableList)
	write.csv (df, outFile, row.names=F, quote=F)

	# Write Prediction Ability file
	outFile = addLabel (fileCrossValPredictions, "GEBVs")
	write.csv (predictionsKFoldsTable, outFile, quote=F, row.names=F)

	return (predictionsKFoldsTable)
}

#-------------------------------------------------------------
# Add label to filename and new extension (optional)
#-------------------------------------------------------------
addLabel <- function (filename, label, newExt=NULL)  {
	nameext = strsplit (filename, split="[.]")
	name    = nameext [[1]][1] 
	if (is.null (newExt))
		ext     = nameext [[1]][2] 
	else
		ext     = newExt
	newName = paste0 (nameext [[1]][1], "-", label, ".", ext )
	return (newName)
}
#-------------------------------------------------------------
#-------------------------------------------------------------
main <- function () {
	library (dplyr)
	library (ggplot2)
	args = commandArgs (trailingOnly=T)
	resultsDir = args[1]
	gs_outputs (resultsDir)
}

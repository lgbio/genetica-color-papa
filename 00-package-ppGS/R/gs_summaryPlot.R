#' Summary plot for cross validation.
#' 
#' Creates a summary ploy by comparing predictive abilities for each phenotype and each model
#'
#' @param predictionsKFoldsFile Filename with data from cross validation: GEBVs for each fold and each repetition.
#' @return None. Writes two files into working dir: data table (.csv) and boxplot graphics (.pdf).
#' @import ggplot2
#' @import dplyr 
#'
#' @export
#-------------------------------------------------------------
#--- Create output tables and boxplots for best model for each trait
#-------------------------------------------------------------
gs_summaryPlot <- function (predictionsKFoldsFile) {
	predictionsKFoldsTable = read.csv (predictionsKFoldsFile)

	message (">>> Creating output plots and tables for best model ...")

	# Summary plot boxplot best models for traits
	meansTable = predictionsKFoldsTable %>% group_by (Traits, Models) %>% 
		summarize (Means=mean(Predictive_ability), .groups="keep") 

	# Get best model for each trait
	bestTable = meansTable %>% group_by (Traits) %>% slice_max (Means) %>% slice_head

	# Create table with correlations from best models
	getBestCorr <- function (trait, model) 
		return (filter (predictionsKFoldsTable, Traits==trait, Models==model))
	bestCorrs = do.call (rbind, mapply (getBestCorr, bestTable$Traits, bestTable$Models, SIMPLIFY=F))

	traitColors = unlist (lapply(1:8, function(x) rep(x,3))) 
	nTraits     = length (unique (bestCorrs$Traits))
	traitColors = traitColors [1: nTraits] 
	outPlotname   = "out-CrossValidation-ModelsTraits.pdf"
	ggplot (bestCorrs, aes(x=HCL_component, y=Predictive_ability)) + ylim (0,1) +
		    geom_boxplot (alpha=0.3, fill=traitColors) + 
			theme(axis.text.x=element_text(angle=0, hjust=1)) + labs (title="Best GS models for color traits in tetraploid potato") + 
			stat_summary(geom='text', label=round (bestTable$Means,2) , fun=max, vjust = -1, size=3, col="blue") +
			stat_summary(geom='text', label=bestTable$Models, fun=min, vjust = 1.5, angle=0,size=3, col="blue") + 
			facet_wrap (~Prefix, ncol=8) 
	ggsave (outPlotname, width=11)
	write.csv (bestCorrs, "out-CrossValidation-ModelsTraits.csv", row.names=F)
}


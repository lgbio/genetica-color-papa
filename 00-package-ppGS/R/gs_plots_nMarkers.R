#-------------------------------------------------------------
# Function for ploting predictive ability varying number of markers
#-------------------------------------------------------------
main <- function () {
	library (ppGS)
	library (ggplot2)
	library (dplyr)
	args = c ("outs-varying-nMarkers")
	gs_plots_varying_nmarkers (args [1])
}
#-------------------------------------------------------------
#-------------------------------------------------------------
#' Create summaries when the number of markers varies. 
#'
#' Create one summary table and two plots: one per base trait and one per hcl trait.
#' @param inputDir  Dir name with results from different GS runs.
#' @return None. Writes table (.csv) and plots (.pdfs) to files.
#' @import dplyr ggplot2
#' @export
gs_plots_varying_nmarkers <- function (inputDir) {
	#library (dplyr)
	#library (ggplot2)

	# Create summary table
	dirList = list.dirs (inputDir, recursive=F)
	summaryTableList = list()
	for (dir in dirList) {
		print (dir)
		n = as.numeric (gsub (paste0(inputDir,"/out"), "", dir))
		hclData = read.csv (sprintf ("%s/%s", dir, "out-CrossValidation-kfolds-GEBVs.csv"))
		summaryTable = data.frame (nMarkers=n, createSummaryTableKFolds (hclData))
		summaryTableList = append (summaryTableList, list (summaryTable))
	}
	summaryTableAll = do.call ("rbind", summaryTableList)
	summaryTableAll
	write.csv (summaryTableAll, "out-GS-crossval-varying-nmarkers.csv", row.names=F)

	# Create plot per base trait
	ggplot (summaryTableAll, aes(x=nMarkers, y=Means, group=HCLs)) + #ylim (0,1) +
			geom_line (aes (color=HCLs)) +
			geom_point (aes (color=HCLs)) +
			labs (title="Genomic selection varying the number of markers per base trait") +
			ylab ("Predictive ability") + xlab ("Number of markers") +
			#scale_x_discrete (name = "Number of Markers", limits = c ("5", "15", "25", "35", "50", "100")) +
			scale_x_continuous(breaks = c(5,15,25,35,50,100), labels = c("5","15","25","35","50","100"),
				limits = c(0,100)) +
			facet_wrap (~Traits, ncol=4) 

	ggsave ("out-GS-crossval-varying-nmarkers-basetraits.pdf", width=11)

	# Create plot per hcl trait
	ggplot (summaryTableAll, aes(x=nMarkers, y=Means, group=Traits)) + 
			#ylim (0,1) +
			geom_line (aes (color=Traits)) +
			geom_point (aes (color=Traits)) +
			labs (title="Genomic selection varying the number of markers per HCL trait") +
			ylab ("Predictive ability") + xlab ("Number of markers") +
			scale_x_continuous(breaks = c(5,15,25,35,50,100), labels = c("5","15","25","35","50","100"),
				limits = c(0,100)) +
			facet_wrap (~HCLs, ncol=3) 


	ggsave ("out-GS-crossval-varying-nmarkers-hcltraits.pdf", width=10, height=7)
}

#-------------------------------------------------------------
#-------------------------------------------------------------
#inputDir = "outs"
#library (dplyr)
#library (ggplot2)
#gs_plots_varying_nmarkers (inputDir)


#-------------------------------------------------------------
#--- Create output tables and boxplots for best model for each trait
#--- Return a summary table with hcl trait, baseTrait, hclComp, best model, and means
#-------------------------------------------------------------
createSummaryTableKFolds <- function (HCLTable) {
	message (">>> Creating summary table of kfolds: hclTrait, baseTrait, hclComp, model, means ...")

	# Summary plot boxplot best models for traits
	meansTable = HCLTable %>% dplyr::group_by (Traits, Models) %>% 
		summarize (Means=mean(get("Predictive_ability")), .groups="keep") 

	# Get best model for each trait
	bestTable = meansTable %>% group_by (Traits) %>% slice_max (Means) %>% slice_head

	# Create table with correlations from best models
	getBestCorr <- function (trait, model) 
		return (dplyr::filter (HCLTable, Traits==trait, Models==model))
	bestCorrs = do.call (rbind, mapply (getBestCorr, bestTable$Traits, bestTable$Models, SIMPLIFY=F))

	# Add base trait and HLC componente to summary table
	baseTraitList = c ()
	hclCompList   = c ()
	for (hclTrait in bestTable$Traits){
		baseTraitList = c (baseTraitList, strsplit (hclTrait, "[.]")[[1]][1])
		hclCompList   = c (hclCompList, strsplit (hclTrait, "[.]")[[1]][3])
	}
	summaryTable = data.frame (bestTable[,1], baseTraitList, hclCompList, bestTable[,2], bestTable [,3])
	names (summaryTable) = c("HCLTraits", "Traits", "HCLs", "Models", "Means")
	write.csv (summaryTable, "out-CrossValidation-GEBVs-Summary.csv", row.names=F)
	return (summaryTable)
}

#-------------------------------------------------------------
#-------------------------------------------------------------
#main ()

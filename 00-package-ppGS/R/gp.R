#' Create summaries when the number of markers varies. 
#'
#' Create one summary table and two plots: one per base trait and one per hcl trait.
#' @param inputDir  Dir name with results from different GS runs.
#' @return None. Writes table (.csv) and plots (.pdfs) to files.
#' @import dplyr ggplot2
#' @export
gs_plots_varying_markers <- function (inputDir) {
	#library (dplyr)
	#library (ggplot2)

	# Create summary table
	dirList = list.dirs (inputDir, recursive=F)
	summaryTableList = list()
	for (dir in dirList) {
		print (dir)
		n = as.numeric (gsub ("outs/out", "", dir))
		hclData = read.csv (sprintf ("%s/%s", dir, "out-CrossValidation-kfolds-GEBVs.csv"))
		summaryTable = data.frame (nMarkers=n, createSummaryTableKFolds (hclData))
		summaryTableList = append (summaryTableList, list (summaryTable))
	}
	summaryTableAll = do.call ("rbind", summaryTableList)
	summaryTableAll
	write.csv (summaryTableAll, "out-kfolds-summary.csv", row.names=F)

	# Create plot per base trait
	ggplot (summaryTableAll, aes(x=nMarkers, y=Means, group=HCLs)) + ylim (0,1) +
			geom_line (aes (color=HCLs)) +
			geom_point (aes (color=HCLs)) +
			labs (title="Genomic selection varying the number of markers per base trait") +
			ylab ("Predictive ability") + xlab ("Number of markers") +
			facet_wrap (~Traits, ncol=4) 
			#theme(axis.text.x=element_text(angle=0, hjust=1)) + labs (title="Best GS models for color traits in tetraploid potato") + 
			#stat_summary(geom='text', label=round (bestTable$Means,2) , fun=max, vjust = -1, size=3, col="blue") +
			#stat_summary(geom='text', label=bestTable$Models, fun=min, vjust = 1.5, angle=0,size=3, col="blue") + 
			#facet_wrap (~Prefix, ncol=8) 

	ggsave ("out-kfolds-summary-traits.pdf", width=11)

	# Create plot per hcl trait
	ggplot (summaryTableAll, aes(x=nMarkers, y=Means, group=Traits)) + ylim (0,1) +
			geom_line (aes (color=Traits)) +
			geom_point (aes (color=Traits)) +
			labs (title="Genomic selection varying the number of markers per HCL trait") +
			ylab ("Predictive ability") + xlab ("Number of markers") +
			facet_wrap (~HCLs, ncol=3) 

			#theme(axis.text.x=element_text(angle=0, hjust=1)) + labs (title="Best GS models for color traits in tetraploid potato") + 
			#stat_summary(geom='text', label=round (bestTable$Means,2) , fun=max, vjust = -1, size=3, col="blue") +
			#stat_summary(geom='text', label=bestTable$Models, fun=min, vjust = 1.5, angle=0,size=3, col="blue") + 
			#facet_wrap (~Prefix, ncol=8) 

	ggsave ("out-kfolds-summary-hcls.pdf", width=11)
}



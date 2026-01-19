#' Create a comparison plot for potato color phenotypes.
#'
#' Function that creates boxplots comparing HCL color components phenotypes (e.g. Heritability or Prediction_ability).
#'
#' @param HCLTable Table with information from cross validation: GEBVs for each fold and each repetition.
#' @param outputDir Output dir for results.
#' @param columnName Column name with predictive ability values (Y axis).
#' @param COLUMNS List with table column names for "phenos", "models", and "prefix"
#' @param TITLES  List with X, Y, and MAIN titles
#' @return None. Plot is saved to a PDF file
#' @import ggplot2
#' @import dplyr 
#' @export
#gs_plotsHCL <- function (HCLFilename, outputDir, X_TITLE, Y_TITLE, TITLE) {
gs_plotsHCL <- function (HCLFilename, outputDir, COLUMNS, TITLES) {

	HCLTable = read.csv (HCLFilename)
	nBaseTraits  = length (unique (HCLTable [,ncol(HCLTable)]))
	message (">>> Creating output plots and tables for best model ...")

    TraitColumn      = COLUMNS$trait
    ModelColumn      = COLUMNS$model
    PrefixColumn     = COLUMNS$prefix
    PredictionColumn = COLUMNS$prediction
    ComponentColumn  = COLUMNS$component

	# Summary plot boxplot best models for traits
	meansTable = HCLTable %>% 
        dplyr::group_by (across (all_of (c(TraitColumn, ModelColumn)))) %>% 
		summarize (Means=mean(.data[[PredictionColumn]], .groups="keep")) 
		#summarize (Means=mean(.data[[TITLES$Y]], .groups="keep")) 

	# Get best model for each trait
	bestTable = meansTable %>% 
        group_by (across (all_of (c(TraitColumn)))) %>% slice_max (Means) %>% slice_head

    view (bestTable)

	# Create table with correlations from best models
	getBestCorr <- function (trait, model) 
		return (dplyr::filter (HCLTable, .data [[TraitColumn]]==trait, .data [[ModelColumn]]==model))

	bestCorrs = do.call (rbind, mapply (getBestCorr, bestTable[[TraitColumn]], bestTable[[ModelColumn]], SIMPLIFY=F))

	traitColors = unlist (lapply(1:nBaseTraits, function(x) rep(x,3))) 
	nTraits     = length (unique (bestCorrs[[TraitColumn]]))
	traitColors = traitColors [1: nTraits] 

	# Scale label text sizes
	
	`%||%` <- function(a, b) if (!is.null(a)) a else b

	current_theme <- theme_get()
	scaled_theme <- theme(
	  plot.title   = element_text(size = (current_theme$plot.title$size %||% 11) * 1.2),
	  axis.title.x = element_text(size = (current_theme$axis.title.x$size %||% 10) * 1.2),
	  axis.title.y = element_text(size = (current_theme$axis.title.y$size %||% 10) * 1.2),
	  axis.text.x  = element_text(size = (current_theme$axis.text.x$size %||% 9) * 1.2),
	  axis.text.y  = element_text(size = (current_theme$axis.text.y$size %||% 9) * 1.2),
	  strip.text   = element_text(size = (current_theme$strip.text$size %||% 9) * 1.2)
	)	
	
    view (bestCorrs)
	#ggplot (bestCorrs, aes(x=get(TITLES$X), y=get(TITLES$Y))) + ylim (0,1) +
	ggplot (bestCorrs, aes(x=HCL_component, y=Predictive_ability)) + ylim (0,1) +
		    geom_boxplot (alpha=0.3, fill=c(1)) + 

		    theme_bw() +	# This changes to a white background with light gray gridlines
			theme (axis.text.x=element_text(angle=0, hjust=1),
                panel.background = element_rect(fill = "white"),  # Ensures pure white background
                plot.background = element_rect(fill = "white"),   # White background for entire plot
                panel.grid.major = element_line(color = "grey90"), # Lighter grid lines
                panel.grid.minor = element_line(color = "grey95")  # Even lighter minor grid lines
            ) + 
			scaled_theme +
			labs (title=TITLES$MAIN, x=TITLES$X, y=TITLES$Y) + 
			stat_summary(geom='text', label=round (bestTable$Means,2) , fun=max, vjust = -1, size=3, col="blue") +
			#stat_summary(geom='text', label=bestTable$Models, fun=min, vjust = 1.5, angle=0,size=3, col="blue") + 
			facet_wrap (~ .data[[PrefixColumn]], ncol=nBaseTraits) 

    #outPlotname    = file.path (outputDir, gsub (".csv", "-PLOT.pdf", basename (HCLFilename)))
	ggsave ("outPlotname.pdf", height=5, width=11)
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

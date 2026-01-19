#' Create files for single genomic selection
#'
#' This function creates both parameters file and training/target datasets for genotype/phenotypes 
#' to use them in genomic selection with single phenotype (gs_single).
#' @param outputDir Directory name where files will be saved.
#' @param nGenos (Optional) Number of markers to include in genotype dataset. If empty, use all markers.
#' @param nPhenos (Optional) Number of phenotypes to include in genotype dataset. If empty, use all phenotypes.
#' @return None
#' @export
gs_testsingle <- function (outputDir, nGenos=NA, nPhenos=NA) {
	createDir (outputDir)
	setwd (outputDir)
	files = gs_data (sprintf ("%s_datasources", outputDir))
	dataDir = sprintf ("%s_datasets", outputDir)
	gs_split (files$geno, files$phenos, dataDir, nGenos, nPhenos=1)
	gs_configfile (type="single", dataDir)
}

#----------------------------------------------------------
# Utility for create dir, if it exists, it is renamed old-XXX
#----------------------------------------------------------
createDir <- function (newDir) {
	checkOldDir <- function (newDir) {
		if (dir.exists (newDir) == T) {
			oldDir = sprintf ("%s/old-%s", dirname (newDir), basename (newDir))
			if (dir.exists (oldDir) == T) checkOldDir (oldDir)
			file.rename (newDir, oldDir)
		}
	}
	checkOldDir (newDir)
	dir.create (sprintf (newDir))
}


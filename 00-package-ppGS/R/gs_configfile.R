#' Creates config file 
#' Creates a template for a config file to be used in genomic selection with multiple models and multiple phenotypes.
#'
#' @param type "multi" or "single" for multiple phenotype or single phenotype configuration file. 
#' @param configDir (Optional) Directory where filenames are locates. Default dir is "."
#' @return None 
#' @export
gs_configfile <- function (type, configDir=".") {
	if (type=="multi") {
		filename="config-multi.yml"
		sink (filename)
			cat (sprintf ("trainingGenoFile  : %s/training-geno.csv  # Obligatory file\n", configDir))
			cat (sprintf ("trainingPhenoFile : %s/training-pheno.csv # Obligatory file\n", configDir))
			cat (sprintf ("targetGenoFile    : %s/target-geno.csv    # Obligatory file\n", configDir))
			cat (sprintf ("targetPhenoFile   : %s/target-pheno.csv   # Optional file \n", configDir))
			cat (sprintf ("markersFile       : NULL                  # Optional file with markers (SNP column) to select for using in GS\n", configDir))
			cat ("nMarkers : 50    # Number of markers to select from markers file\n")
			cat ("methods  : BA BB # Models to test, the complete models are : BA BB BC BL BRR GBLUP EGBLUP RKHS\n")
			cat ("nTimes   : 1     # Number of times for cross validation, use 5 or more for better results (take long time)\n")
			cat ("nFolds   : 2     # Number of folds for cross validation, use 5 or more for better results (take long time)\n") 
			cat ("nCores   : 4     # Number of processing cores for GS, it dependes of the machine where it is running\n")
		sink()
	}else {
		filename="config-single.yml"
		sink (filename)
			cat (sprintf ("trainingGenoFile  : %s/training-geno.csv  # Obligatory file\n", configDir))
			cat (sprintf ("trainingPhenoFile : %s/training-pheno.csv # Obligatory file\n", configDir))
			cat (sprintf ("targetGenoFile    : %s/target-geno.csv    # Obligatory file\n", configDir))
			cat (sprintf ("targetPhenoFile   : %s/target-pheno.csv   # Optional file \n", configDir))
			cat ("phenoName : PriFlor.LCH.H         # Name of the phenotype \n")
			cat ("methods   : BA BB # Models to test, the complete models are : BA BB BC BL BRR GBLUP EGBLUP RKHS\n")
			cat ("nTimes    : 1     # Number of times for cross validation, use 5 or more for better results (take long time)\n")
			cat ("nFolds    : 2     # Number of folds for cross validation, use 5 or more for better results (take long time)\n") 
			cat ("nCores    : 4     # Number of processing cores for GS, it dependes of the machine where it is running\n")
		sink()
	}
}

#------------------------------
# Main only for test in command line
#------------------------------
main <- function () {
	args = commandArgs (trailingOnly=T)
	gs_configFile (args[1])
}
#main ()

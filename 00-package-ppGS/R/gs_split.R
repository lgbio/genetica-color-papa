#' Split genotype/phenotype datasets into training and target datasets
#'
#' Function that splits a complete genotype/phenotypes data into training/target datasets.
#'
#' @param genoFile   Filename of genotype data matrix (rows for markers, columns for samples)
#' @param phenosFile Filename of fenotype data table (rows for samples, columns for phenotypes)
#' @param outputDir  The directory where the results will be saved.
#' @param nGenos     (Optional) Number of markers to use from genotype data matrix. If empty all markers are used. 
#' @param nPhenos    (Optional) Number of phenotype to use from phenotype data table. If empty all phenotypes are used. 
#' 
#' @return None, results are saved to the output directory
#'
#' @export
#--------------------------------------------------------
# Splits dataset into training and testing datasets\n
#--------------------------------------------------------
gs_split <- function (genoFile, phenosFile, outputDir, nGenos=NA, nPhenos=NA) {
	set.seed (1)
	message (">>> Splitting genotype into training/target datasets...")
	
	createDir (outputDir)
	message ("Current dir: ", getwd())

	# Preprocess geno and feno 
	geno  = preprocessingGenotype (genoFile, outputDir, nGenos)
	pheno = preprocessingPhenos (phenosFile, outputDir, nPhenos, imputed=T)	

	# Select only common individuals from geno and pheno
	genoIndividuals        = colnames (geno)
	phenoIndividuals       = rownames (pheno)
	commonIndividuals      = intersect (genoIndividuals, phenoIndividuals) 
	uncommonIndividuals    = setdiff (genoIndividuals, commonIndividuals)

	# Select samples for training and testing (20% for testing)
	nTesting               = round (0.2*length (commonIndividuals))
	targetIndividuals      = sample (commonIndividuals, size=nTesting)
	trainingIndividuals    = setdiff (commonIndividuals, targetIndividuals)

	#

	# Create datasets
	trainingGeno           = cbind (Samples=trainingIndividuals, t (geno  [,trainingIndividuals]))
	trainingPhenos         = cbind (Samples=trainingIndividuals, pheno [trainingIndividuals,,drop=F])
	targetGeno             = cbind (Samples=targetIndividuals, t (geno  [,targetIndividuals]))
	targetPhenos           = cbind (Samples=targetIndividuals, pheno [targetIndividuals,,drop=F])
	rownames (targetPhenos) = targetIndividuals

	# Write datasets
	write.csv (trainingGeno,   sprintf ("%s/%s", outputDir, "training-geno.csv"), quote=F, row.names=F)
	write.csv (trainingPhenos, sprintf ("%s/%s", outputDir, "training-pheno.csv"), quote=F, row.names=F)
	write.csv (targetGeno,     sprintf ("%s/%s", outputDir, "target-geno.csv"), quote=F, row.names=F)
	write.csv (targetPhenos,   sprintf ("%s/%s", outputDir, "target-pheno.csv"), quote=F, row.names=F)
}

#--------------------------------------------------------------------
# Preprocesses genotype: Currently, only get a subset (nGenos) of genos
# No imputation with mice package as it takes too long
#--------------------------------------------------------------------
preprocessingGenotype <- function (genoFile, outputDir, nGenos=NA) {
	geno = read.csv (genoFile, check.names=F, row.names=1)
	#out = mice (geno, m=5, maxit=0, meth="pmm", seed=500)
	#genoImputed = complete (out, 1)
	if (is.na (nGenos)){
		nGenos = nrow (geno)-1
	}
	newGeno = geno [1:nGenos,]
	return (newGeno)
}

#--------------------------------------------------------------------
# Preprocesses phenotype by imputing and removing outliers
#--------------------------------------------------------------------
preprocessingPhenos <- function (phenosFilename, outputDir, nPhenos=NA, imputed=T) {
	#---------------- local functions---------------#
	#---- Impute null or NA values using "mice" package
	imputePhenos <- function (df) {
		view (df)
		init  = mice::mice(df, maxit=0) 
		meth  = init$method
		predM = init$predictorMatrix
		set.seed(103)
		imputed     = mice::mice(df, method=meth, predictorMatrix=predM, m=5, print=FALSE)
		impComplete = mice::complete (imputed)
		return (impComplete)
	}

	#---- Drop 1% of extreme values
	dropOutliers <- function (df, percent=0.01) {
		n          = nrow (df)
		nDrop      = floor (percent*n)
		dfOutliers = df [nDrop:(n-nDrop),]
		return (dfOutliers)
	}
	#---------------- local functions---------------#
	data         = read.csv (phenosFilename, check.names=T)
	if (is.na (nPhenos))
		nPhenos = ncol (data)-1

	names        = data [,1, drop=F]
	phenos       = data [,-1, drop=F]
	phenos       = phenos [,1:nPhenos, drop=F]

	if (imputed==FALSE) { 
		phenos = read.csv (phenosFilename, check.names=F, row.names=1)
		return (phenos[,1:nPhenos, drop=F])
	}

	dataImputed  = imputePhenos (phenos)
	view (dataImputed)
	dataImputed  = data.frame (names, dataImputed)
	phenos = dropOutliers (dataImputed)
	rownames (phenos) = phenos [,1]
	phenos [,1] = NULL

	phenosFilename = basename (addLabel (phenosFilename, "IMPUTED"))
	write.csv (phenos, sprintf ("%s/%s", outputDir, phenosFilename), row.names=T)
	return (phenos)
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

#----------------------------------------------------------
# Create dir, if it exists the it is renamed old-XXX
#----------------------------------------------------------
createDir <- function (newDir) {
	checkOldDir <- function (newDir) {
		name  = basename (newDir)
		path  = dirname  (newDir)
		if (dir.exists (newDir) == T) {
			oldDir = sprintf ("%s/old-%s", path, name)
			if (dir.exists (oldDir) == T) 
				checkOldDir (oldDir)
			file.rename (newDir, oldDir)
		}
	}
	checkOldDir (newDir)
	dir.create (newDir)
}
#-------------------------------------------------------------
# Main for test from command line
#-------------------------------------------------------------
main <- function () {
	#library (mice)
	USAGE="USAGE: gs_split.R <genotype matrix filename> <phenotype table filename> <output dir> [number of genotypes] [number of phenotypes]"
	args  = commandArgs(trailingOnly = TRUE)
	gs_split (args[1], args[2], args[3], args[4], args[5])
}

#main ()

utils::globalVariables(c("Traits", "Models", "y", "Predictive_ability", "SNP", "Means", "HCL_component"))
#' Genomic selection for multiple phenotypes and multiple models
#'
#' Runs in parallel genomic selection for multiple phenotypes and multiple models.
#' 
#' @param configFile Filename of the parameter file (YAML Format).  It includes:
#'        trainingGenoFile: Filename of training genotype data,
#'        trainingPhenoFile: Filename of training phenotype data, 
#'        trainingGenoFile: Filename of target genotype table, 
#'        targetPhenoFile: Filename of target phenotype table (opt.),
#'        markersFile: Filename of selected markers to test from genotype (opt.),
#'        nMarkers: Number of selected markers to test,
#'        methods: Methods or models to test in GS,
#'        nTimes: Number of independent replicates for the cross-validation, 
#'        nFolds: Number of folds for the cross-validation,
#'        nCores: Number of cores to use.
#'
#' @param outputDir The directory where the results will be saved.
#'
#' @return None
#'
#' @import parallel
#' @import yaml
#' @import ggplot2
#' @import dplyr 
#'
#' @export 
gs_multi <- function (configFile, outputDir) {
    message ("+++ Running gs_multi...")
	#-------------------------------------------------------------
	# Multiple phenotype / Multiple models genomic selection
	#-------------------------------------------------------------
	params     = loadConfigurationParameters (configFile)
	phenosDir  = "phenotypes"

	# Read/Processing input files and set directories
	init       = initFilesAndWorkingDir (outputDir, params,  configFile)
	files      = createInputFilesGS (init$data, init$files, params, phenosDir) 

	# Create parameters file for each phenotype to run GS 
	paramsPhenoFiles = createParamFileForEachPhenotype (files, params, phenosDir)

	# Run in parallel GS
	dummy=mclapply (paramsPhenoFiles, run_gs_single, init$outputDir, mc.cores=params$nCores)
	message ("LOG Errors: ")
	print (dummy)

	# Organize result files
	organizeOutputFiles (init$outputDir)

	# Create summary plots and tables
	message ("Creating predictionsKFoldsTable...")
	resultsDir              = "results"
	predictionsKFoldsTableFilename  = createSummaryTablesGS (resultsDir)
    return (predictionsKFoldsTableFilename)
	#createSummaryPlot (predictionsKFoldsTable)
}
#-------------------------------------------------------------
# Wrapped function for gs_single
#-------------------------------------------------------------
run_gs_single <- function (configFile, outputDir) {
	setwd (outputDir)
	out = gs_single (configFile)
	return (out)
}

#-------------------------------------------------------------
# Create parallelGS parameter file for each phenotype
#-------------------------------------------------------------
createParamFileForEachPhenotype <- function (files, params, phenosDir) {

	paramsPhenoFiles = c()
	for (i in 1:length(files$phenoNames)) {
		phenoName     = files$phenoNames [i] 
		trainingGeno  = files$trainingGeno [i]
		targetGeno    = files$targetGeno [i]
		trainingPheno = files$trainingPhenos [i]
		targetPheno   = files$targetPhenos [i]

		configFile = sprintf ("%s/%s-params.yml", phenosDir, phenoName)
		sink (configFile)
		cat (sprintf ("trainingGenoFile  : %s\n", trainingGeno))
		cat (sprintf ("trainingPhenoFile : %s\n", trainingPheno))
		cat (sprintf ("targetGenoFile    : %s\n", targetGeno))
		cat (sprintf ("targetPhenoFile   : %s\n", targetPheno))
		cat (sprintf ("phenoName         : %s\n", phenoName))
		cat (sprintf ("useGwasMarkers    : %s\n", params$useGwasMarkers))    # Use only best markers from GWAS results
		cat (sprintf ("gwasResultsPath   : %s\n", params$gwasResultsPath))
		cat (sprintf ("methods           : %s\n", params$methods))
		cat (sprintf ("nTimes            : %s\n", params$nTimes))
		cat (sprintf ("nFolds            : %s\n", params$nFolds))
		cat (sprintf ("nCores            : %s\n", params$nCores))
		sink ()

		paramsPhenoFiles = c (configFile, paramsPhenoFiles)
	}
	return (paramsPhenoFiles)
}

#-------------------------------------------------------------
# Move files to directories: data, phenotypes, and results
#-------------------------------------------------------------
organizeOutputFiles <- function (outputDir) {
	setwd (outputDir)
	dir.create ("data")

	for (f in list.files (pattern=".csv"))
		file.rename (f, sprintf ("data/%s", f))
	for (f in list.files (pattern=".yml"))
		file.rename (f, sprintf ("data/%s", f))
	for (f in list.files (pattern=".log"))
		file.rename (f, sprintf ("data/%s", f))

	dataFiles = setdiff (list.files (), c("data", "phenotypes"))
	dir.create ("results")
	for (f in dataFiles) 
		file.rename (f, sprintf ("results/%s", f))
}

#-------------------------------------------------------------
# Read input files and load them to "data" var
# Create output directory and copy into it genotype and phenotype data
#-------------------------------------------------------------
initFilesAndWorkingDir <- function (outputDir, params, configFile, createNewDir=T) {
	data = c()
	if (!all (sapply(c(params$trainingGenoFile, params$trainingPhenoFile,
					  params$targetGenoFile), file.exists))) {
		message ("Any of the input files do not exist!")
		quit ()
	}

	data$trainingGeno   = read.csv (params$trainingGenoFile)
	data$trainingPhenos = read.csv (params$trainingPhenoFile)
	data$targetGeno     = read.csv (params$targetGenoFile)
	data$nMarkers       = params$nMarkers
	# Check if gwas markers was given
	if (is.null(params$markersFile)) 
		data$gwasMarkers = NULL
	else
		# Get N random markers, only for testing 
		if (params$markersFile=="RANDOM") {
			randomMarkers     = unlist (sample (colnames (data$trainingGeno[,-1]), data$nMarkers))
			data$trainingGeno = cbind (SAMPLES=data$trainingGeno[,1], data$trainingGeno [,randomMarkers])
			data$targetGeno   = cbind (SAMPLES=data$targetGeno[,1], data$targetGeno [,randomMarkers])
		}else
			data$gwasMarkers = read.csv (params$markersFile)

	# Check if target pheno was given
	if (is.null (params$targetPhenoFile)) 
		data$targetPhenos = NULL
	else
		data$targetPhenos = read.csv (params$targetPhenoFile)

	files = c()
	files$trainingGenoFile    = addLabel (basename (params$trainingGenoFile), "TRAINING")
	files$targetGenoFile      = addLabel (basename (params$targetGenoFile), "TARGET")
	files$trainingPhenoFile   = addLabel (basename (params$trainingPhenoFile), "TRAINING")

	# Create and set global working dir where all outputs will be saved
	if (createNewDir)
		createDir (outputDir)
	file.copy (configFile, outputDir)
	setwd (outputDir)
	outputDir = getwd()  # Full path

	write.csv (data$trainingGeno,   files$trainingGenoFile, quote=F, row.names=F)
	write.csv (data$targetGeno,     files$targetGenoFile, quote=F, row.names=F)
	write.csv (data$trainingPhenos, files$trainingPhenoFile, quote=F, row.names=F)
	if (!is.null (params$targetPhenoFile)) {
		files$targetPhenoFile  = addLabel (basename (params$targetPhenoFile), "TARGET")
		write.csv (data$targetPhenos, files$targetPhenoFile, quote=F, row.names=F)
	}

	return (list (data=data, files=files, outputDir=outputDir))
}

#-------------------------------------------------------------
# Create genotype and phenotype files for each phenotype
#-------------------------------------------------------------
createInputFilesGS <- function (data, files, params, phenosDir) {
	# Write one file by trait
	phenoNames            = colnames (data$trainingPhenos)[-1];#view (phenoNames)
	trainingGenoFileList  = c()
	targetGenoFileList    = c()
	trainingPhenoFileList = c()
	targetPhenoFileList   = c()
	createDir (phenosDir)
	for (i in 1:length(phenoNames)) {
		message (phenoNames [i])
		phenoName = phenoNames [i] 

		# Get all or selected (GWAS) markers
		genoFiles        = getMarkersForGS (data$gwasMarkers, data$nMarkers, data$trainingGeno, data$targetGeno, files$trainingGenoFile, files$targetGenoFile,  phenoName, phenosDir)

		trainingGenoFile = genoFiles$training
		targetGenoFile   = genoFiles$target

		# Get single phenotype
		trainingPheno           = data$trainingPhenos [, c(1,1+i)] # ["Samples", "XXX"]
		trainingPhenoFilename   = sprintf ("%s/%s-trainingPheno.csv", phenosDir, phenoName)
		write.csv (trainingPheno, trainingPhenoFilename, quote=F, row.names=F)

		trainingGenoFileList    = c (trainingGenoFileList, trainingGenoFile)
		targetGenoFileList      = c (targetGenoFileList, targetGenoFile)
		trainingPhenoFileList   = c (trainingPhenoFileList, trainingPhenoFilename)

		targetPhenoFile = "NULL"
		if (!is.null (params$targetPhenoFile)) {
			targetPheno         = data$targetPhenos [, c(1,1+i)]   # ["Samples", "XXX"]
			targetPhenoFile     = sprintf ("%s/%s-targetPheno.csv", phenosDir, phenoName)
			write.csv (targetPheno, targetPhenoFile, quote=F, row.names=F)
		}
		targetPhenoFileList = c(targetPhenoFileList, targetPhenoFile)
	}

	outList = list(trainingGenoFileList, targetGenoFileList, trainingPhenoFileList, files$targetGenoFile, targetPhenoFileList, phenoNames)
	names (outList) = c("trainingGeno", "targetGeno", "trainingPhenos", "targetGeno", "targetPhenos", "phenoNames")
	return (outList)
}

#-------------------------------------------------------------
# Use ALL or GWAS markers 
# Return the genotype filename with ALL or GWAS markers
# gwasMarkers is a file with N best markers for trait sorted by agroScore
#-------------------------------------------------------------
getMarkersForGS <- function (gwasMarkers, nMarkers, trainingGeno, targetGeno, trainingGenoFile, targetGenoFile, phenoName, phenosDir) {
	if (is.null (gwasMarkers)) {
		message (">>> GS with ALL markers...")
	}else {
		message (">>> GS with GWAS markers for: ", phenoName )
		bestGwasMarkers   = gwasMarkers$SNP [gwasMarkers$TRAIT==phenoName]
		nMarkersRatio     = floor (nMarkers * length (bestGwasMarkers) / 100)
		message ("+++ Ratio of ", nMarkers, "% markers for ", phenoName, ": ", nMarkersRatio)
		bestGwasMarkers   = head (gwasMarkers$SNP [gwasMarkers$TRAIT==phenoName], nMarkersRatio)
		view (bestGwasMarkers)
		#bestGwasMarkers   = unique (gwasMarkers$SNP [gwasMarkers$TRAIT==phenoName])
		#bestGwasMarkers   = getNScoredBestSNPs (gwasMarkers, phenoName, nMarkers)

		allMarkers        = colnames (trainingGeno);
		commonMarkers  = intersect (allMarkers, bestGwasMarkers) 
		#view (allCommonMarkers)
        #nCommonMarkers    = allCommonMarkers [1:nMarkers]
        #nCommonMarkers    = allCommonMarkers [1:nMarkers]
		#view (nCommonMarkers)
        

		trainingGeno      = data.frame (SAMPLES=trainingGeno [,1], trainingGeno [,commonMarkers])
		trainingGenoFile  = sprintf ("%s/%s-trainingGeno.csv", phenosDir, phenoName)
		write.csv (trainingGeno, trainingGenoFile, quote=F, row.names=F)

		targetGeno      = data.frame (SAMPLES=targetGeno [,1], targetGeno [,commonMarkers])
		targetGenoFile  = sprintf ("%s/%s-targetGeno.csv", phenosDir, phenoName)
		write.csv (targetGeno, targetGenoFile, quote=F, row.names=F)

		#message ("...END...");quit()
	}
	return (list (training=trainingGenoFile, target=targetGenoFile))
}

#-----------------------------------------------------------
# Obsolote as the gwasMarkers
# Sort best SNPs from multiple tools using a own score
#-----------------------------------------------------------
old_getNScoredBestSNPs <- function (gwasMarkers, phenoName, nMarkers) {
	scoresTable = gwasMarkers [gwasMarkers[,1]==phenoName,]

	# Add Count of SNPs between groups
	dfCountSNPs  = data.frame (add_count (scoresTable, SNP, sort=F, name="scoreShared")); 

	# Score GC
	scoreSign   = ifelse (scoresTable$SIGNIFICANCE, 1, 0)
	scoreGC     = 1 - abs (1-scoresTable$GC)
	scoreSNPs   = data.frame (add_count (scoresTable, SNP))$n
	SCORE_SNP   = 0.5*scoreGC + 0.3*scoreSNPs + 0.2*scoreSign

	dataSNPsScores     = select (data.frame (scoresTable, SCORE_SNP), SNP, SCORE_SNP, everything())
	sortedGwasMarkers  = dataSNPsScores [order (-SCORE_SNP),] %>% distinct (SNP, .keep_all=T)
	outBestFilename    = sprintf ("outBestTrait-%s.csv", phenoName)
	write.csv (sortedGwasMarkers, outBestFilename, row.names=F)

	nSNPsScores      = nrow (sortedGwasMarkers)

	if (nMarkers > nSNPsScores)
		nMarkers = nSNPsScores
	#bestGwasMarkers    = sort (levels (scoresTable$SNP))
	bestGwasMarkers  = unlist (sortedGwasMarkers [1:nMarkers,1])
	return (bestGwasMarkers)
}

#----------------------------------------------------------
# Read configuration parameters
#----------------------------------------------------------
loadConfigurationParameters <- function (configFile) {
	params = yaml.load_file (configFile)
	return (params)
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
	system (sprintf ("mkdir %s", newDir))
}

#-------------------------------------------------------------
# Create summary tables from GS results of each phenotype
# Create a GEBVs table and kfolds table, the last is returned to be graphed
#-------------------------------------------------------------
createSummaryTablesGS <- function (inputDir) {
	message ("++++ Creating table from cross validations predictive abilities...")

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

    message ("++++++++++++++++++ Returning filename:", outFile)
	return (outFile)
}

#-------------------------------------------------------------
#--- Create output tables and boxplots for best model for each trait
#-------------------------------------------------------------
createSummaryPlot <- function (predictionsKFoldsTable) {
	nBaseTraits = length (unique (HCLTable [,1]))
	message (">>> Creating output plots and tables for best model ...")

	# Summary plot boxplot best models for traits
	meansTable = predictionsKFoldsTable %>% dplyr::group_by (Traits, Models) %>% 
		summarize (Means=mean(Predictive_ability), .groups="keep") 

	# Get best model for each trait
	bestTable = meansTable %>% group_by (Traits) %>% slice_max (Means) %>% slice_head

	# Create table with correlations from best models
	getBestCorr <- function (trait, model) 
		return (dplyr::filter (predictionsKFoldsTable, Traits==trait, Models==model))
	bestCorrs = do.call (rbind, mapply (getBestCorr, bestTable$Traits, bestTable$Models, SIMPLIFY=F))

	traitColors = unlist (lapply(1:nBaseTraits, function(x) rep(x,3))) 
	nTraits     = length (unique (bestCorrs$Traits))
	traitColors = traitColors [1: nTraits] 
	outPlotname   = "out-CrossValidation-ModelsTraits.pdf"
	ggplot (bestCorrs, aes(x=HCL_component, y=Predictive_ability)) + ylim (0,1) +
		    geom_boxplot (alpha=0.3, fill=traitColors) + 
			theme(axis.text.x=element_text(angle=0, hjust=1)) + 
			labs (title="Best GS models for LCH components in potato color traits") +
			stat_summary(geom='text', label=round (bestTable$Means,2) , fun=max, vjust = -1, size=3, col="blue") +
			stat_summary(geom='text', label=bestTable$Models, fun=min, vjust = 1.5, angle=0,size=3, col="blue") + 
			facet_wrap (~Prefix, ncol=nBaseTraits) 
	ggsave (outPlotname, width=11)
}

#-------------------------------------------------------------
#-------------------------------------------------------------

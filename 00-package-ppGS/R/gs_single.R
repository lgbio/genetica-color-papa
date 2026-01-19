utils::globalVariables (c("EN", "brewer.pal"))

main <- function () {
	library (ppGS)
	library (yaml)
	library (parallel)
	gs_single ("SecBrote.C-params.yml", "out")
}

#' Genomic selection for single phenotype and multiple models 
#'
#' Function that performs genomic selection by evaluating in parallel multiple models on a single phenotype.
#' Models are evaluated through cross-validation, then the best model is selected, 
#' and GEBVs are calculated using this model. 
#' Results are generated in both tables and graphics, and they include: 
#' GEBVs for target and training datasets, graphic comparison of predicted ability for each model 
#' (correlation between actual and predicted values), and comparison of execution times between models.
#'
#' @param paramsFile Filename of the parameter file (YAML Format).
#' @param outputDir  (Optional) Directory name where results will be saved. If empty, phenotype name is used. 
#'
#' @return TRUE if no errors or FALSE otherwise.
#'
#' @import utils
#' @import grDevices
#' @import graphics
#' @importFrom stats complete.cases
#' @import parallel
#' @import yaml
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
gs_single <- function (paramsFile, outputDir="") {
	message (">>>>>>>>>>>>>>>>>>>> Running GS for parameters file: ", paramsFile, "<<<<<<<<<<<<<<<<<<<<<<<")
	#------------------------------------------------------
	# Single phenotype / Multiple models genomic selection
	#------------------------------------------------------

	# load input files with or without all markers
	out    = readGenoPhenoTrainingTargetData (paramsFile)
	data   = out$data
	params = out$params

	errorOutput = TRUE
	tryCatch ({
		trainingGeno  = data$trainingGeno
		trainingPheno = data$trainingPheno
		targetGeno    = data$targetGeno
		targetPheno   = data$targetPheno

		# Create and change to output dir
		if (outputDir=="")
			outputDir = data$phenoName
		createDir (outputDir)
		setwd (outputDir)

		res = runCrossValidation (data, params)
		listCVs = res$listCVs
		params$methods = res$methods
		createOutputsCrossValidation (listCVs, params$phenoName)
		bestMethod = getBestGenomicPredictionMethod (listCVs, params$methods)

		predictGEBVsTraining (trainingGeno, trainingPheno, trainingGeno, res$methods, params$phenoName, params$nCores)
		predictGEBVsTarget (trainingGeno, trainingPheno, trainingGeno, trainingPheno, bestMethod, params$phenoName, "TRAINING")
		predictGEBVsTarget (trainingGeno, trainingPheno, targetGeno, targetPheno, bestMethod, params$phenoName, "TARGET")

		flush(stderr())
		#BglrMethods = c("BA", "BB", "BC", "BL","BRR","RKHS")
		#METHODS = c("GBLUP","EGBLUP","RR","LASSO","EN","BRR","BL","BA", "BB", "BC","RKHS","RF","SVM")
		#AllMethods = c("GBLUP","EGBLUP","RR","LASSO","EN","BRR","BL","BA", "BB", "BC","RKHS","RF","SVM","BRNN")
	}, error = function (e){
		message (sprintf ("Error in gs_single: In trait %s ", outputDir), e)
		errorOutput <<- FALSE
	})
	return (errorOutput)
}

#-----------------------------------------------------------------------------------------------------------
# Compare predicted value with diffetent methods
#-----------------------------------------------------------------------------------------------------------
predictGEBVsTarget <- function (trainingGeno, trainingPheno, targetGeno, targetPheno,  bestMethod, TRAIT, PREFIX) {
	message (">>> Predicting GEBVs for best model using target data...")

	# Run predictions target
	params             = list ("gtrn"=trainingGeno, "ptrn"=trainingPheno, "gtrg"=targetGeno)
	predictedGEBVs     = predictGEBVsFunction (bestMethod, params)

	if (is.null (targetPheno)) 
		predictedGEBVsTable =  data.frame (predictedGEBVs, row.names = rownames (targetGeno))
	else {
		targetPheno         = as.data.frame (targetPheno)
		noNAs               = row.names (targetPheno [complete.cases (targetPheno),,drop=F])
		predicted           = data.frame (noNAs, PREDICTED=predictedGEBVs [noNAs], row.names=1)
		predictedGEBVsList  = list (predicted, ACTUAL=targetPheno [noNAs,])

		# Write outputs to table and histograms
		mainText            = sprintf ("%s GEBVs for '%s' trait, using '%s' model", PREFIX, TRAIT, bestMethod)
		outFile             = sprintf ("out-GEBVs-BestModel-%s", PREFIX)
		predictedGEBVsTable = as.data.frame (predictedGEBVsList)
		rownames (predictedGEBVsTable) = noNAs
		message (">>> Creating outputs for Actual vs. Predicted GEBVs...", TRAIT)
		graphicPredictedOutputs (predictedGEBVsList, mainText, outFile, predictedGEBVsTable)
	}
	colnames (predictedGEBVsTable)[1] = sprintf ("GEBVs(%s)", bestMethod)
	predictedGEBVsTable = cbind (SAMPLES=rownames (predictedGEBVsTable), predictedGEBVsTable)
	write.csv (predictedGEBVsTable, "out-GEBVs-BestModel-TARGET.csv", row.names=F,quote=F)
}

#-----------------------------------------------------------------------------------------------------------
# Compare predicted value with diffetent methods
#-----------------------------------------------------------------------------------------------------------
predictGEBVsTraining <- function (trainingGeno, trainingPheno, targetGeno, METHODS, TRAIT, NCORES) {
	message (">>> Predicting GEBVs for different models only using training data...")
	params = list ("gtrn"=trainingGeno, "ptrn"=trainingPheno, "gtrg"=targetGeno)

	# Run predictions in parallel
	message (">>> METHODS:\n", METHODS)
	predictedGEBVsList = mclapply (METHODS, predictGEBVsFunction, params, mc.cores=NCORES, mc.silent=T)

	# Write outputs to table and histograms
	outFile             = "out-GEBVs-AllModels-TRAINING"
	mainText            = sprintf ("Comparison of GEBVs from different models for '%s' trait", TRAIT)
	predictedGEBVsTable = as.data.frame (predictedGEBVsList, col.names=METHODS)
	message (">>> Creating outputs for Predicted GEBVs from different methods...")
	graphicPredictedOutputs (predictedGEBVsList, mainText, outFile, predictedGEBVsTable)
}

#----------------------------------------------------------
# Read geno/pheno trainig/target data
# Also, use selected or all markers, if "gwasResultsPath" is set or not.
#----------------------------------------------------------
readGenoPhenoTrainingTargetData <- function (paramsFile) {
	message (">>> Reding geno/pheno datasets...")
	# Load parameters file and make a list from methods parameter
	params         = yaml.load_file (paramsFile)
	params$methods = strsplit (params$methods, "[ ]+")[[1]]

	# Read arguments and load data
	trainingGeno     = read.csv (params$trainingGenoFile, check.names=F, row.names=1)
	trainingPhenoTbl = read.csv (params$trainingPhenoFile, check.names=F)
	targetGeno       = read.csv (params$targetGenoFile, check.names=F, row.names=1)
	phenoName        = params$phenoName
	# Check if target pheno was given
	if (is.null (params$targetPhenoFile)) 
		targetPheno = NULL
	else {
		targetPhenoTbl      = read.csv (params$targetPhenoFile, check.names=F)
	  	targetPheno         = targetPhenoTbl [,2]
	  	names (targetPheno) = targetPhenoTbl [,1]
	}

	# Convert phenotype table to vector
	trainingPheno         = trainingPhenoTbl [,2]
	names (trainingPheno) = trainingPhenoTbl [,1]
	trait                 = colnames (trainingPhenoTbl)[2]

	data = list (trainingGeno, trainingPheno, targetGeno, targetPheno, trait, phenoName)
	names (data) = c ("trainingGeno", "trainingPheno", "targetGeno", "targetPheno", "trait", "phenoName")

	return (list (params=params, data=data))
}
#-------------------------------------------------------------
# Wrapper for BWGS function "predict" that predics GEBVs
#-------------------------------------------------------------
predictGEBVsFunction <- function (method, params) {
	message (">>>>>>>>>>", method)
	predictedGEBVs  = bwgs.predict (predict.method=method,
	                                geno_train=params$gtrn,pheno_train=params$ptrn, geno_target=params$gtrg,
	                                MAXNA=0.2,MAF=0.05,geno.reduct.method="NULL",reduct.size="NULL",
	                                r2="NULL",pval="NULL",MAP="NULL",geno.impute.method="MNI")
	df = predictedGEBVs [,1]
	names (df) = rownames (params$gtrg)
	return (df)
}

#-------------------------------------------------------------
# Prediction ability comparison of methods
# a.Execute prediction methods
#-------------------------------------------------------------
runCrossValidation <- function (data, params) {
	trainingGeno  = data$trainingGeno
	trainingPheno = data$trainingPheno
	METHODS       = params$methods
	nTimes        = params$nTimes
	nFolds        = params$nFolds
	NCORES        = params$nCores
	TRAIT         = params$phenoName

	message ("\t>>> Running cross validation using nTimes=",nTimes, " nFolds=", nFolds)
	message ("\t>>> Methods=", METHODS)

	crossValidationFunction <- function (method) {
		message (">>> ", method)

        # Create out dir for trait + method
		createDir (method)
		setwd (method)

		cvOutput = NA
		tryCatch ({
			cvOutput <- bwgs.cv (trainingGeno, trainingPheno, predict.method= method,
	                       geno.impute.method="mni", nFolds=nFolds, nTimes=nTimes )
	  		return (cvOutput)
	    }, error=function (e) {
			message (sprintf ("Error running crossValidationFunction: method: %s ", method), e)
			view (trainingGeno)
			view (trainingPheno)
			return (NA)
		})
	}

	listCVs         = mclapply (METHODS, crossValidationFunction, mc.cores=NCORES, mc.silent=T)
	names (listCVs) = METHODS
	# Remove NAs items from listCVs and METHOS
	METHODS = METHODS [!is.na (listCVs)]
	listCVs = listCVs [!is.na (listCVs)]
	return (list (listCVs=listCVs, methods=METHODS))
}

#-------------------------------------------------------------
# Create output tables and graphics from Cross Validation results
#-------------------------------------------------------------
createOutputsCrossValidation <-function (listCVs, TRAIT) {
	message ("\t>>> Get info from cross validation...")
	predictedAbilitiesList   = list()
	predictedTimesList       = list()
	predictedAbilitiesKFolds = list()

	for (cvObject in listCVs) {
		predictedAbilitiesList   = append (predictedAbilitiesList, list (cvObject$cv))
		predictedTimesList       = append (predictedTimesList, list (cvObject$summary[[3]]))
		predictedAbilitiesKFolds = append (predictedAbilitiesKFolds, list (cvObject$summary$KFoldsTable))
	}
	predAbilTable            = as.data.frame (predictedAbilitiesList)
	predTimesTable           = as.data.frame (predictedTimesList)
	predAbilKFoldsTable      = do.call ("rbind", predictedAbilitiesKFolds)

	#predictedAbilitiesMeans = mclapply (listCVs, function (cvObject ){return (cvObject$summary[[1]])})
	#bestMethod          = METHODS [which.max (predictedAbilitiesMeans)]

	# Write outputs: tables and boxplots
	nModels   = ncol (predAbilTable)
	allColors = unlist (lapply (c ("Set1", "Set2", "Set3"), function (x) brewer.pal (n=8, x)))
	COLORS    = allColors [1:nModels]

	message ("\t>>> Outputs for predicted ability by methods...")
	outFilename = "out-CrossValidation"

	#write.csv (predAbilTable, paste0 (outFilename, "-means.csv"), quote=F, row.names=F)
	#pdf (paste0 (outFilename, "-means.pdf"), width=13, height=13)
	#mainText = paste0 ("Predictive ability by models\nTrait: ", TRAIT)
	#boxplot(predAbilTable, xlab="Prediction model",  ylab="Predictive ability", main=mainText, col=COLORS)
	#dev.off()

	message ("\t>>> Outputs for times for each method...")
	write.csv (predTimesTable, paste0 (outFilename, "-times.csv"), quote=F, row.names=F)
	pdf (paste0 (outFilename, "-times.pdf"), width=13, height=13)
	mainText = paste0 ("Average time by models\nTrait: ", TRAIT)
	barplot (unlist (predictedTimesList), main=mainText, ylab="Time (mins)", xlab="Methods", col=COLORS)
	dev.off()

	# Create kfolds table with PA for each fold and time
	message ("\t>>> Outputs for predicted ability by methods correlations by fold...")
	write.csv (predAbilKFoldsTable, paste0 (outFilename, "-kfolds.csv"), quote=F, row.names=F)
	pdf (paste0 (outFilename, "-kfolds.pdf"), width=13, height=13)
	mainText = paste0 ("Cross validation by models \nTrait: ", TRAIT)
	boxplot (Predictive_ability~Models, data=predAbilKFoldsTable, main=mainText, col=COLORS)
	dev.off()
}
#-----------------------------------------------------------------------------------------------------------
# c. Compare predicted value varying training size
#-----------------------------------------------------------------------------------------------------------
comparePredictionAbilitiesBySize <- function (listCVs, trainingGeno, trainingPheno, METHODS, nTimes, nFolds, TRAIT, NCORES, NCHUNKS) {
	message (">>> Predicting values by size")
	nSamples  = length (trainingPheno)
	N         = NCHUNKS   # number of chunks
	chunkSize = floor (nSamples / N)
	sizes     = seq (chunkSize, nSamples, chunkSize)
	sizes     = sizes [1:N]

	message ("\t>>> Get info from cross validation")
	predictedAbilitiesMeans = mclapply (listCVs, function (cvObject ){return (cvObject$summary[[1]])}, mc.silent=T)
	bestMethod          = METHODS [which.max (predictedAbilitiesMeans)]

	#--- Local Function -------------------------------------------------------------------
	cvFunSize <- function (sampleSize, method) {
	  cvOutput <-bwgs.cv (sample.pop.size=sampleSize,
	                      predict.method=method,
	                      geno=trainingGeno, pheno=trainingPheno, pop.reduct.method="RANDOM",
	                      geno.impute.method="mni", nFolds=nFolds, nTimes=nTimes )
	  return (cvOutput$cv)
	} #-------------------------------------------------------------------------------------

	listCVsSize = mclapply (sizes[1:(N-1)], cvFunSize, bestMethod, mc.cores=NCORES, mc.silent=T)
	listCVsSize = append (listCVsSize, list(listCVs[[bestMethod]]$cv))
	names (listCVsSize) = paste0("N=", sizes)

	message ("\t>>> Writing outputs for predicted value by size...")
	outFilename = "out-predicted-value-by-size-bestMethod"
	CompareSize = as.data.frame (listCVsSize, col.names=paste0("N=",as.character(sizes)),check.names=F)
	write.csv (CompareSize, paste0 (outFilename, ".csv"), quote=F, row.names=F)

	pdf (paste0 (outFilename, ".pdf"), width=13, height=13)
	mainText = paste0 ("Effect of training population size\n Trait: ", TRAIT, ", Method: ", bestMethod )
	boxplot(CompareSize,xlab="Training POP size",ylab="Predictive avility", main=mainText)
	dev.off()
}

#-------------------------------------------------------------
# Graphic histograms and scatterplots for GEBVx
#-------------------------------------------------------------
graphicPredictedOutputs <- function (predictedGEBVsList, mainText, outFilename, predictedGEBVsTable, PREFIX="") {
	## put histograms on the diagonal
	panel.hist <- function(x, ...) {
	  usr = par("usr"); on.exit(par(usr))
	  par(usr = c(usr[1:2], 0, 1.5) )
	  h      = hist(x, plot = FALSE)
	  breaks = h$breaks; nB <- length(breaks)
	  y      = h$counts; y <- y/max(y)
	  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
	}
	## put (absolute) correlations on the upper panels, with size proportional to the correlations.
	panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
	  usr = par("usr"); on.exit(par(usr))
	  par(usr = c(0, 1, 0, 1))
	  r   = abs(cor(x, y))
	  #txt = format(c(r, 0.123456789), digits = digits)[1]
	  txt = format(r, digits = digits)
	  txt = paste0(prefix, txt)
	  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
	  text(0.5, 0.5, txt, cex = cex.cor * r)
	}

	write.csv (predictedGEBVsTable, paste0 (outFilename, ".csv"), quote=F, row.names=T)
	pdf (paste0 (outFilename, ".pdf"), width=13, height=13)
	pairs(predictedGEBVsTable,lower.panel = panel.smooth,
	      upper.panel = panel.cor, diag.panel=panel.hist, main=mainText)
	dev.off()
}

#-------------------------------------------------------------
# Get Best Method
#-------------------------------------------------------------
getBestGenomicPredictionMethod <- function (listCVs, METHODS) {
	predictedAbilitiesMeans = mclapply (listCVs, function (cvObject ){return (cvObject$summary[[1]])}, mc.silent=T)
	bestMethod          = METHODS [which.max (predictedAbilitiesMeans)]
	return (bestMethod)
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

#----------------------------------------------------------
#----------------------------------------------------------
#main ()

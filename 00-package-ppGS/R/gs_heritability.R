#' Calculates and plots heritabilities for HCL penotypes.
#' 
#' @param genoFile   Filename of genotype matrix.
#' @param phenosFile Filename of phenotypes table.
#' @param outputDir  Dirname where outputs will be saved.
#' @return None. Generate two files: a table (CSV) with heritabilities for each trait
#'         and a summary plot (PDF) comparing heritabilities from all traits.
#' @import parallel
#' @import dplyr
#' @import ggplot2
#' @import BGLR
#' @export
gs_heritability <- function (genoFile, phenosFile, outputDir) {
	data       = readProcessGenoPheno (genoFile, phenosFile)
	geno       = data$geno
	pheno      = data$pheno 
	traitNames = colnames (pheno)

	message (">>> Calculating heritabilities...")
	createDir (outputDir)
	createDir (sprintf ("%s/tmp", outputDir))

	h2Results = mclapply (traitNames, calculateHeritabilityTrait, 
						  geno, pheno, outputDir, mc.cores=4)
	hclFile   = createHCLTable (traitNames, h2Results, outputDir)
	gs_plotsHCL (hclFile, outputDir, "Heritability", 
				 "Comparison of mean heritabilities for phenotypes")
}

#-------------------------------------------------------------
# Plot heritability by multi boxplots
#-------------------------------------------------------------
createHCLTable <- function (traitNames, h2Results, outputDir) {
	message (">>> Creating HCL table...")

	h2Table = data.frame ()
	for (i in 1:length (h2Results)) {
		Prefix        = strsplit (traitNames[i], "[.]")[[1]][1]
		HCL_component = strsplit (traitNames[i], "[.]")[[1]][2]
		Models        = "BA"
		Traits        = traitNames [i]
		Heritability  = h2Results [[i]]
		tmpDF  = data.frame (Prefix, HCL_component, Models, 
							 Traits, Heritability,stringsAsFactors=F)
		h2Table = rbind (h2Table, tmpDF)
	}
	hclTable = h2Table [order (h2Table$Traits),]
	hclFile  = sprintf ("%s/out-%s-HCL-Comparison-TABLE.csv", outputDir, "Heritability")
	write.csv (hclTable, hclFile , quote=F, row.names=F)
	return (hclFile)
}

#-------------------------------------------------------------
# Calculata heritability for a single trait using BayesA method
# It can be improved by testing all methods and selecting best ones for trait
#-------------------------------------------------------------
calculateHeritabilityTrait <- function (traitName, geno, pheno, outputDir) {
	message (">>> Heritability for ", traitName)
	# Get Trait and Geno for BGLR
	X = scale(geno)/sqrt(ncol(geno))
	Y = pheno [, traitName]
	 
	# Estimation using variance components
	MODEL  = "BayesA"
	PREFIX = paste0 (outputDir, "/tmp/", MODEL, "-", traitName,"-")
	#varUFilename = paste0 (PREFIX, "ETA_1_varB.dat")
	varUFilename = paste0 (PREFIX, "ETA_1_ScaleBayesA.dat")
	varEFilename = paste0 (PREFIX, "varE.dat")

	h2_0 = .5
	eta  = list(list(X=X,model=MODEL,saveEffects=T))
	fm   = BGLR(y=Y,ETA=eta, nIter=5000, burnIn=1000, verbose=F, saveAt=PREFIX)
	varU = scan(varUFilename)
	varE = scan(varEFilename)
	h2   = varU/(varU+varE)

	return (h2)
}

#-------------------------------------------------------------
# Read genotype and phenotype, transform, impute, and filter by MAF
#-------------------------------------------------------------
readProcessGenoPheno <- function (genotypeFile, phenotypeFile) {
	# Read data, get geno matrix, and process it for BGRL
	genoCCC  = read.csv (genotypeFile, check.names=F, row.names=1);view(genoCCC)
	phenoCCC = read.csv (phenotypeFile, check.names=F,  row.names=1);view (phenoCCC)
 
	genoCCC  = t (genoCCC);view (genoCCC)
   
	samplesGeno   = rownames (genoCCC);#view (samplesGeno)
	samplesPheno  = rownames (phenoCCC);#view (samplesPheno)
	samplesCommon = intersect (samplesGeno, samplesPheno); #view (samplesCommon)
	view (samplesCommon)

	geno  = genoCCC  [samplesCommon,]; #view (geno)
	pheno = phenoCCC [samplesCommon,]; #view (pheno)

	view (geno)
	view (pheno)

	#-------------------------------------------------------------
	# Impute NA alleles
	#-------------------------------------------------------------
	imputeGenotype <- function (M) {
		message(">>> Missing marker data imputed with population mode...")
		impute.mode <- function(x) {
			ix <- which(is.na(x))
			if (length(ix)>0) 
				x[ix] <- as.integer(names(which.max(table(x))))
			return(x)
		}
		missing <- which(is.na(M))
		if (length(missing)>0) {
			M <- apply(M,2,impute.mode)
		}
		return (M)
	}
	genoImputed = imputeGenotype (geno);#view (genoImputed)

	#-------------------------------------------------------------}
	# MAF genotype
	#-------------------------------------------------------------}
	MAFGenotype <- function (M, thresholdMAF) {
		message (">>> Checking minor allele frecuency, MAF=", thresholdMAF)
		calcMAF <- function(x) {
			ploidy=4
			AF <- mean(x, na.rm=T) / ploidy;
			MAF <- ifelse(AF > 0.5,1-AF,AF)
		}
		MAF <- apply(M, 2, calcMAF)
		polymorphic <- which(MAF>thresholdMAF)
		M <- M[,polymorphic]
	}
	genoMAF = MAFGenotype (genoImputed, 0.1);#view(genoMAF)

	return (list(geno=genoMAF, pheno=pheno))
}
		
#-------------------------------------------------------------
# Plot heritabilities without grouping HCL components
#-------------------------------------------------------------
plotHeritabilities <- function (traitNames, h2Results) {
	message (">>> Plotting heritabilities...")

	if (is.null (h2Results))
		h2Table = read.csv ("heritability-values.csv")
	else {
		h2Table = data.frame ()
		for (i in 1:length (h2Results)) 
			h2Table = rbind (h2Table, data.frame (Heritability=h2Results[[i]], Traits=traitNames [i]))	
		write.csv (h2Table, "heritability-values.csv", quote=F, row.names=F)
	}

	# Write means by trait
	h2Means=group_by(h2Table, Traits) %>% summarize(Means=mean (Heritability), .groups="drop") %>% arrange (-Means)
	write.csv (h2Means, "heritability-means-traits.csv", quote=F, row.names=F)


	n <- length (h2Results)
	library (RColorBrewer)
	COLORS=brewer.pal(n = 8, name = "Set2")
	palette = c()
	for (c in COLORS) 
		palette = c(palette, rep (c,3))
	#palette <- sample (colors(distinct=T), n/3) 

	pdf ("heritability-traitsColors-boxplots.pdf", width=15)
		par(mar=c(6, 4.1, 4.1, 2.1))
		bx=boxplot (Heritability~Traits, h2Table, col=palette, xaxt="n", xlab="", 
				 main="Heritability for color traits in potato")

		# Write medians at top 
		bxMedians = round (bx$stats[3,], 2) 
		text(x =  seq_along(traitNames), y = par("usr")[4]+0.02, srt = 0, adj = 0.5, labels = bxMedians, xpd = TRUE, col="red")

		axis(1, at=1:length(traitNames), labels = FALSE)
		labels = gsub ("LCH-","", traitNames)
		text(x =  seq_along(labels), y = par("usr")[3]-0.03, srt = 45, adj = 1, labels = labels, xpd = TRUE)
		# Limit lines
		abline(h=0.3,lty=2,col="red",lwd=2)
		text (x=0.0, y=0.28,labels=c("Low (0.0~0.3)"), adj=0, col="red")
		text (x=0.0, y=0.58,labels=c("Medium (0.31~0.6)"), adj=0, col="blue")
		abline(h=0.6,lty=2,col="blue",lwd=2)
		text (x=0.0, y=0.98,labels=c("High (0.61~1.0)"), adj=0, col="darkgreen")
		abline(h=1.0,lty=2,col="darkgreen",lwd=2)

		mtext ("Traits", side=1, col="black", line=5)
	dev.off()


	h2_0 = 0.5
	pdf ("heritability-scatterplots.pdf", width=11)
	op=par (mfrow = c(4,6))
		for (i in 1:length(h2Results)) {
			traitName = traitNames [i]
			h2 = h2Results [[i]]
		 	plot(h2,type='o',cex=.5,col=4, main=traitName);
			abline(h=c(h2_0,mean(h2[-c(1:200)])),lty=2,col=c(1,2),lwd=2)
		}
	par (op)
	dev.off()

	pdf ("heritability-histograms.pdf", width=11)
	op=par (mfrow = c(4,6))
		for (i in 1:length(h2Results)) {
			traitName = traitNames [i]
			h2 = h2Results [[i]]
	 		hist (h2, main=traitName)
		}
	par (op)
	dev.off()

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
#-------------------------------------------------------------
# Main only for testing
#-------------------------------------------------------------
main <-  function() {
	#library (ppGS)
	source ("lglib14.R")
	source ("gs_plotsHCL.R")
	library (parallel)
	library (dplyr)
	library (ggplot2)
	library (BGLR)
	args = c ("inputs/genos.csv", "inputs/phenos.csv", "inouts")
	
	gs_heritability (args [1], args [2], args [3])
	gs_plotsHCL ("inouts/out-Heritability-HCL-Comparison-TABLE.csv", "inouts", "Heritability", 
				 "Comparison of mean heritabilities for phenotypes")
}

#main ()


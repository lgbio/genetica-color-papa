#!/usr/bin/env Rscript

# --- tell renv the project BEFORE it loads ---
Sys.setenv(RENV_PROJECT = "/home/lg/BIO/agrosavia/genetica-color-papa")

if (!requireNamespace("renv", quietly = TRUE)) stop("renv not available")
#renv::activate(project = "/home/lg/BIO/agrosavia/genetica-color-papa")
renv::load ()

cat("R version:", R.version.string, "\n")
cat("Library paths:\n")
print(.libPaths())

source ("lglib14.R")

library (reshape2) # melt
library (ggplot2)

# Plot densities and outliers of HCL components using density and boxplot plots
# Outliers removed according to Interquartile Range (IQR) method
# Plot densities and outliers after removing outliers

#INPUTFILE = "inputs/out-fenotiposColor-ComponentesLCH.csv"
INPUTFILE = "inputs/fenotiposColor-ComponentesLCH.csv"
createDir ("outputs")

#----------------------------------------------------------
#----------------------------------------------------------
main <- function () {
	source ("lglib14.R")
	args = commandArgs (trailingOnly=T)
	args = c(INPUTFILE)
	
	phenosFile  = args [1]
	plotDensitiesHCLTraits (phenosFile)
	plotOutliersHCLTraits (phenosFile)

	# Remove outliers
	noOutlieresPhenoFile = removeOutliers (phenosFile)
	plotDensitiesHCLTraits (noOutlieresPhenoFile)
	plotOutliersHCLTraits (noOutlieresPhenoFile)
}

#-------------------------------------------------------------
#-------------------------------------------------------------
# Remove outliers according to Interquartile Range (IQR) method
#-------------------------------------------------------------
removeOutliers <- function (phenosFile) {
	phenos    = read.csv (phenosFile);  ;#view (phenos)
	traitNames = names (phenos [,-1]) # remove "Register" column

	rmOutliers <- function (values) {
		outliers = boxplot (values, plot=F)$out
		positions = which (values %in% outliers)
		return (positions)
	}

	for (trait in traitNames) {
		message ("Trait: ", trait)
		values = phenos [,trait]
		positions = rmOutliers (values)
		if (length (positions) > 0)
			message ("Outliers detected: in trait: ", trait, ", length: ", length(positions))
			phenos [positions,trait] = NA
	}
	
	outFile = gsub (".csv", "-noOUTLIERS.csv", phenosFile)
	outFile  = gsub ("inputs/", "outputs/", outFile)
	write.csv (phenos, outFile, row.names = F)
	return (outFile)
}

#----------------------------------------------------------
#----------------------------------------------------------
plotOutliersHCLTraits <- function (phenosFile) {
	message ("plotOutliersHCLTraits...")
	phenos = read.csv (phenosFile)
	hclFile = createHCLTable (phenosFile)

	hclTable    = read.csv (hclFile)
	traitColors = unlist (lapply(1:9, function(x) rep(x,3))) 
	nTraits     = length (unique (hclTable$HCLTrait))
	traitColors = traitColors [1: nTraits] 

	ggplot (hclTable, aes(x=HCL, y=Value)) + #ylim (0,1) +
		    geom_boxplot (alpha=0.5, fill=traitColors) + 
			theme(axis.text.x=element_text(angle=0, hjust=1)) + 
			theme(text = element_text(size = 12)) +
			labs (title="Outliers HCL components", x="HCL Component", y="Value") + 
			facet_wrap (~TRAIT, ncol=nTraits) 

	outFile  = gsub ("inputs/", "outputs/", phenosFile)
	outFile  = gsub (".csv", "-BOXPLOTS.pdf", outFile)
	ggsave (outFile, width=7)
	return (outFile)
}

#-------------------------------------------------------------
#-------------------------------------------------------------
createHCLTable <- function (phenosFile) {
	phenosOrg  = read.csv (phenosFile);  ;#view (phenosOrg)
	phenosRsh  = reshape2::melt (phenosOrg[,-1]);  ;#view (phenosRsh)
	names (phenosRsh) <- c("HCLTrait", "Value")

	trait    = sapply (strsplit (as.character (phenosRsh [,1]),  "[.]"), function (x) x[1])
	hcl      = sapply (strsplit (as.character (phenosRsh [,1]),  "[.]"), function (x) x[2])
	#phenosRsh$HCLTrait = paste0 (trait,".",hcl)

	#hclTrait   = strsplit (as.character (phenosRsh [,1]), "[.]")

	outFile    = gsub (".csv", "-HCLTABLE.csv", phenosFile)
	outFile  = gsub ("inputs/", "outputs/", outFile)
	HCLTable   = cbind (TRAIT=trait, HCL=hcl, phenosRsh)
	write.csv (HCLTable, outFile, row.names =F)
	return (outFile)
}

#----------------------------------------------------------
# Graficar las distribuciones de las variables de los phenos de entrada
#----------------------------------------------------------
plotDensitiesHCLTraits <- function (phenosFile) {
	message ("+++ plotDensitiesHCLTraits...")
	hclFile  = createHCLTable (phenosFile)
	hclTable = read.csv (hclFile)

	# Plot densities
	ggplot (hclTable, aes (x=Value)) +
			geom_density() +
			theme (text = element_text(size = 8), 
				   strip.text = element_text(size = 12)) +
			labs (title="Densities HCL components", x="Values",y="Distribution density ") +
			facet_wrap (~HCLTrait, scale="free", ncol=6) 


	outputFile  = gsub ("inputs/", "outputs/", phenosFile)
	outputFile = gsub (".csv", "-DENSITIES.pdf", outputFile) 
    message ("+++ ", outputFile)
  	ggsave (outputFile, width=10)
}

#----------------------------------------------------------
#----------------------------------------------------------
main ()

suppressPackageStartupMessages(library(synthpop))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(sampling))
suppressPackageStartupMessages(library(partykit))

library (ppGS)
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

mycols <- c("darkmagenta", "turquoise")
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
myseed <- 20190110

# Get datasets
gs_datasets ("sources")

# Change samples names in both geno/pheno
genoFile   = "sources/sources-genotypes-CCC-Andigena-ClusterCall2020.csv"
phenoFile = "sources/sources-phenotypes-CCC-Andigena-ColorTraitsHCL.csv"
geno  = read.csv (genoFile, check.names=F)
pheno = read.csv (phenoFile)

genoSamples = colnames (geno)[-1]
phenoSamples = pheno [,1]

allSamples = union (genoSamples, phenoSamples)
nSamples = length (allSamples)
samples = sapply (1:nSamples, function (x) sprintf ("Sample%.3d", x))
samples.df = data.frame (REG=allSamples, SAMPLE=samples)

for (i in 1:nSamples) {
	org = samples.df [i,1]
	newn = as.character (samples.df [i,2])
	names (geno)[names(geno)==org] = newn
	pheno [pheno[,1]==org,1] = newn
}

names (geno)
pheno[,1]

# Change trait names in pheno
cols = colnames (pheno)[-1]
ncols = c("Samples")
for (i in 1:length (cols)) 
	ncols = c(ncols, sprintf ("Pheno%.2d", i))

print (ncols)
colnames (pheno) = ncols

# Change markers names in geno
markers    = geno [,1]
nMarkers   = length (markers)
newMarkers = sapply (1:nMarkers, function (x) sprintf ("Marker%.4d",x)) 
geno [,1]  = newMarkers

# Write synthetized files
outGeno = addLabel (genoFile, "SYNTH")
write.csv (geno, outGeno, row.names=F)

outPheno = addLabel (phenoFile, "SYNTH")
write.csv (pheno, outPheno, row.names=F)

# Synthesise data
synth.obj <- syn(pheno, seed = myseed)
phenosSyn = synth.obj$syn

# compare the synthetic and original data frames
pdf ("cmp.pdf")
compare(synth.obj, pheno, nrow = 3, ncol = 4, cols = mycols)
dev.off()


#!/usr/bin/env Rscript
#detach("package:ppGS", unload=TRUE)
library (ppGS)
# Get sourde geno and pheno datasets (pheno: componentsLCH)
#gs_datasets ("sources")

# Create configuration file multiple phenotypes
#gs_configfile (type="multi")
set.seed (123)

gs_split ("inouts/genotipo.csv",
		  "inouts/fenotipos.csv",
		  outputDir="data-ADDITIVE")
#, nGenos=1000, nPheno=3)
gs_multi ("config-multi-TEST.yml", "outputs-ADDITIVE")

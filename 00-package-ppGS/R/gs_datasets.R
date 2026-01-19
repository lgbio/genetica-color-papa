#' Write datasets to files
#'
#' Function that writes into CSV files different types of datasets used for genomica selection (training/target/sources phenotypes and genotypes)
#'
#' @param type    String defining the type of dataset to generate: "small", "full", "sources, or "single". 
#'                "small" and "full" create both a configuration file (YAML format) and training/target datasets for genomic selection with multiple phentotypes. The former with a genotype of 1000 markers and a phenotype of 6 traits, and the latest with a genotype of 4640 markers and 24 traits.
#'                "sources" copies original genotype and phenotype datasets used to generate the "small" and "full" datasets.
#'                And "single" creates both a configuration file (YAML format) and training/target datasets for genomic selection with single phenotype.
#' @return None. Files are written into directories with the same name of the type arguments (i.e, "small", "full", "sources", and "single")
#' @export
gs_datasets <- function (type) {
	message (">>> Copying ", type, " datasets to ", getwd())
	if (type=="sources") {
		data ("sourceGeno")
		data ("sourcePhenos")
		createDir (type)
		genoFile   = sprintf ("%s/%s", type, "sources-genotypes-CCC-Andigena-ClusterCall2020.csv")
		write.csv (sourceGeno, genoFile, row.names=F)
		phenosFile = sprintf ("%s/%s", type, "sources-phenotypes-CCC-Andigena-ColorTraitsHCL.csv")
		write.csv (sourcePhenos, phenosFile, row.names=F)
	}else if (type=="small"){
		data ("smallTrGn")
		data ("smallTgGn")
		data ("smallTrPh")
		data ("smallTgPh")
		data ("smallConfig")
		createDir (type)
		message ("small-training-geno.csv")
		write.csv (smallTrGn, sprintf ("%s/%s", type, "small-training-geno.csv"), row.names=F)
		message ("small-target-geno.csv")
		write.csv (smallTgGn, sprintf ("%s/%s", type, "small-target-geno.csv"), row.names=F)
		message ("small-training-pheno.csv")
		write.csv (smallTrPh, sprintf ("%s/%s", type, "small-training-pheno.csv"), row.names=F)
		message ("small-target-pheno.csv")
		write.csv (smallTgPh, sprintf ("%s/%s", type, "small-target-pheno.csv"), row.names=F)
		message ("small-config.yml")
		write.table (smallConfig, sprintf ("%s/%s", type, "small-config.yml"), col.names=F, row.names=F, quote=F)
	}else if (type=="full"){
		data ("fullTrGn")
		data ("fullTgGn")
		data ("fullTrPh")
		data ("fullTgPh")
		createDir (type)
		message ("full-training-geno.csv")
		write.csv (fullTrGn, sprintf ("%s/%s", type, "full-training-geno.csv"), row.names=F)
		message ("full-target-geno.csv")
		write.csv (fullTgGn, sprintf ("%s/%s", type, "full-target-geno.csv"), row.names=F)
		message ("full-training-pheno.csv")
		write.csv (fullTrPh, sprintf ("%s/%s", type, "full-training-pheno.csv"), row.names=F)
		message ("full-target-pheno.csv")
		write.csv (fullTgPh, sprintf ("%s/%s", type, "full-target-pheno.csv"), row.names=F)
		message ("full-config.yml")
		write.table (fullConfig, sprintf ("%s/%s", type, "full-config.yml"), col.names=F, row.names=F, quote=F)
	}else if (type=="single"){
		data ("singleTrGn")
		data ("singleTgGn")
		data ("singleTrPh")
		data ("singleTgPh")
		createDir (type)
		message ("single-training-geno.csv")
		write.csv (singleTrGn, sprintf ("%s/%s", type, "single-training-geno.csv"), row.names=F)
		message ("single-target-geno.csv")
		write.csv (singleTgGn, sprintf ("%s/%s", type, "single-target-geno.csv"), row.names=F)
		message ("single-training-pheno.csv")
		write.csv (singleTrPh, sprintf ("%s/%s", type, "single-training-pheno.csv"), row.names=F)
		message ("single-target-pheno.csv")
		write.csv (singleTgPh, sprintf ("%s/%s", type, "single-target-pheno.csv"), row.names=F)
		message ("single-config.yml")
		write.table (singleConfig, sprintf ("%s/%s", type, "single-config.yml"), col.names=F, row.names=F, quote=F)
	}else
		message ('>>> Error: not valid type. Valid type are: "small", "full", "sources", or "single"')
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

#!/usr/bin/env Rscript

#------------------- For calling from a renv environment -------------------------------------------------
Sys.setenv(RENV_PROJECT = "/home/lg/BIO/agrosavia/genetica-color-papa/14x-SeleccionGenomica-N-Marcadores")
source(file.path(Sys.getenv("RENV_PROJECT"), "renv/activate.R"))
#---------------------------------------------------------------------------------------------------------

# ---- Dependencies ----
library (dplyr)
library (BGLR)
library (parallel)
library (coda)
library (ggplot2)

source ("lglib14.R")
source ("gs_plotsHCL_20.R")
# ---- Defaults (tune as needed) ----
NITER   <- 12000
BURNIN  <- 2000
THIN    <- 10
VERBOSE <- FALSE
MODEL    = "BRR"
DEBUG    = FALSE
NCORES   = 7

main <-  function() {
    rootDir          = getwd ()
    INPUTSDIR        = ifelse (!DEBUG, "inputs/", "inputs/test/")
    phenotypesFile	 = paste0 (INPUTSDIR, "fenotipos.csv")
    genotypeFile	 = paste0 (INPUTSDIR, "genotipo.csv")
    englishNamesFile = paste0 ("inputs/traitnames-BASE.csv")
    
    # Create info table of missingness per trait
    log_trait_missingness (phenotypesFile, "info-trait-missingness.csv")
    quit ()

    #MODELS    = c("RKHS","BRR","BayesA","BayesB","BayesC","BayesRR","BL")
    MODELS    = c("RKHS", "BRR")
    for (model in MODELS) {
        setwd (rootDir)
        OUTPUTDIR = paste0 ("outputs-", model)

        # NOTE: use the loop 'model' here (not "BRR")
        hclTableFilename = gs_heritability (genotypeFile, phenotypesFile, OUTPUTDIR, model, NCORES)

        hclTableFilename_SPANISH = hclTablePreprocessing (hclTableFilename, englishNamesFile, OUTPUTDIR, "SPANISH")
        gs_plotsHCL (hclTableFilename_SPANISH, OUTPUTDIR, "Componentes_CHL",  "Heredabilidad",
                     "Heredabilidad para los 21 componentes de color CHL")

        hclTableFilename_ENGLISH = hclTablePreprocessing (hclTableFilename, englishNamesFile, OUTPUTDIR, "ENGLISH")
        gs_plotsHCL (hclTableFilename_ENGLISH, OUTPUTDIR, "CHL_Components", "Heritability", 
                     "Heritability for the 21 CHL Components")
    }
}

# ===========================================
# Per-trait heritability (posterior chain)
# ===========================================
heritabilityPolyTrait <- function(
  traitName, geno, pheno, ETA, outputDir,
  MODEL = "RKHS",
  nIter = NITER, burnIn = BURNIN, thin = THIN, verbose = VERBOSE
) {
  message("Calculating Heritability for trait: ", traitName, "...")
  y <- pheno[, traitName]
  keep <- !is.na(y)
  y <- y[keep]
  ETAk <- subset_ETA(ETA, keep)

  resultsDir <- file.path(outputDir, traitName)
  createDir(resultsDir)

  fit <- BGLR(
    y = y, ETA = ETAk,
    nIter = nIter, burnIn = burnIn, thin = thin,
    verbose = verbose, saveAt = paste0(resultsDir, "/")
  )

  add_suffix       <- getVarFilenameForModel(MODEL)
  additiveVarFile  <- file.path(resultsDir, paste0("ETA_1_", add_suffix))
  residualVarFile  <- file.path(resultsDir, "varE.dat")

  # No file checks by design—assumes BGLR wrote them
  var_additive <- scan(additiveVarFile)
  var_residual <- scan(residualVarFile)

  h2 <- var_additive / (var_additive + var_residual)
  return (h2)
}

gs_heritability <- function( genoFile, phenosFile, outputDir, MODEL, NCORES=1,
                             min_n_per_trait = 10) {  # set to 0 to disable skipping
  message(">>> \n\nCalculating heritabilities (MODEL = ", MODEL, ")...")
  createDir(outputDir)
  createDir(file.path(outputDir, "tmp"))

  # Preprocess once (align, impute, MAF)
  data  <- readProcessGenoPheno(genoFile, phenosFile)
  geno  <- data$geno
  pheno <- data$pheno

  # 3a) Log missingness for transparency
  log_trait_missingness(pheno, file.path(outputDir, "trait_missingness.csv"))

  # ETA built once and reused for every trait
  ETA <- build_ETA(geno, MODEL)

  # Trait list: numeric phenotype columns
  traitNames_all <- names(pheno)[vapply(pheno, is.numeric, TRUE)]
  if (!length(traitNames_all)) stop("No numeric trait columns found in phenotype table.")

  # 3b) (Optional) filter ultra-sparse traits
  nn <- colSums(!is.na(pheno[, traitNames_all, drop = FALSE]))
  traitNames <- traitNames_all[nn >= min_n_per_trait]
  if (!length(traitNames)) stop("No traits meet min_n_per_trait = ", min_n_per_trait)

  # Compute h2 posteriors per trait
  h2Results <- mclapply(
    traitNames, heritabilityPolyTrait,
    geno = geno, pheno = pheno, ETA = ETA, outputDir = outputDir, MODEL = MODEL,
    mc.cores = NCORES
  )
  names(h2Results) <- traitNames

  # Summaries
  summary_df <- do.call(
    rbind,
    lapply(names(h2Results), function(tr) {
      h2 <- h2Results[[tr]]
      qs <- quantile(h2, c(0.025, 0.5, 0.975), na.rm = TRUE)
      data.frame(
        trait  = tr,
        mean   = mean(h2),
        median = qs[[2]],
        ci2.5  = qs[[1]],
        ci97.5 = qs[[3]],
        n_samp = length(h2),
        model  = MODEL,
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(summary_df) <- NULL

  # Also write a compact summary CSV for quick inspection
  out_csv <- file.path(outputDir, paste0("heritability_summary_", MODEL, ".csv"))
  write.csv(summary_df, out_csv, row.names = FALSE)

  # IMPORTANT: pass the MODEL used here (don’t read the global)
  hclTableFilename = createHCLTable (traitNames, h2Results, outputDir, MODEL)

  return (hclTableFilename)
}

createHCLTable <- function (traitNames, h2Results, outputDir, model_used) {
  message (">>> Creating HCL table...")

  h2Table = data.frame ()
  for (i in seq_along(h2Results)) {
    Prefijo          = strsplit (traitNames[i], "[.]")[[1]][1]
    Componentes_CHL  = strsplit (traitNames[i], "[.]")[[1]][2]
    Modelos          = model_used          # <- use the argument passed in
    Fenotipos        = traitNames [i]
    Heredabilidad    = h2Results [[i]]
    tmpDF  = data.frame (Prefijo, Componentes_CHL, Modelos,
                         Fenotipos, Heredabilidad, stringsAsFactors = FALSE)
    h2Table = rbind (h2Table, tmpDF)
  }
  hclTable = h2Table [order (h2Table$Fenotipos),]
  hclTableFilename  = sprintf ("%s/out-%s-HCL-Comparison-TABLE.csv", outputDir, "Heritability")
  write.csv (hclTable, hclTableFilename , quote = FALSE, row.names = FALSE)
  return (hclTableFilename)
}

log_trait_missingness <- function(pheno, out_csv) {
  nn  <- colSums(!is.na(pheno))
  nax <- colSums(is.na(pheno))
  df  <- data.frame(trait = names(nn), n_non_missing = as.integer(nn), n_missing = as.integer(nax))
  df  <- df[order(df$n_non_missing, decreasing = TRUE),]
  write.csv(df, out_csv, row.names = FALSE)
  message("Wrote trait missingness report to: ", out_csv)
  invisible(df)
}

# ===========================================
# Your preprocessing (kept intact)
# ===========================================
# Read genotype and phenotype, transform, impute, and filter by MAF
readProcessGenoPheno <- function (genotypeFile, phenotypeFile) {
  genoCCC  = read.csv (genotypeFile,  check.names = FALSE, row.names = 1)
  phenoCCC = read.csv (phenotypeFile, check.names = FALSE, row.names = 1)

  genoCCC  = t(genoCCC)

  samplesGeno   = rownames(genoCCC)
  samplesPheno  = rownames(phenoCCC)
  samplesCommon = intersect(samplesGeno, samplesPheno)

  geno  = genoCCC [samplesCommon, ]
  pheno = phenoCCC[samplesCommon, ]

  # Impute NA alleles by population mode (for polyploid dosages)
  imputeGenotype <- function (M) {
    message(">>> Missing marker data imputed with population mode...")
    impute.mode <- function(x) {
      ix <- which(is.na(x))
      if (length(ix) > 0)
        x[ix] <- as.integer(names(which.max(table(x))))
      x
    }
    if (any(is.na(M))) M <- apply(M, 2, impute.mode)
    M
  }
  genoImputed = imputeGenotype(geno)

  # MAF filter for tetraploid dosages (0..4)
  MAFGenotype <- function (M, thresholdMAF) {
    message (">>> Checking minor allele frequency, MAF = ", thresholdMAF)
    calcMAF <- function(x) {
      ploidy <- 4
      AF  <- mean(x, na.rm = TRUE) / ploidy
      ifelse(AF > 0.5, 1 - AF, AF)
    }
    MAF <- apply(M, 2, calcMAF)
    polymorphic <- which(MAF > thresholdMAF)
    M[, polymorphic, drop = FALSE]
  }
  genoMAF = MAFGenotype(genoImputed, 0.10)

  return (list(geno = genoMAF, pheno = pheno))
}

#--------------------------------------------------------------------------------
# Change to english names, remove two traits
# Original trait names in spanish. Titles in english
#--------------------------------------------------------------------------------
hclTablePreprocessing <- function (hclTableFilename, englishNamesFile, outputDir, outType) {
	mappings <- read.csv (englishNamesFile) # Load the name mappings
  
	df <- read.csv (hclTableFilename)       # Load the target table

  	# Remove specific traits. 
	`%notin%` <- Negate(`%in%`)
	df <- df %>% filter(Prefijo %notin% c('CBaya', 'CPulpa'))

	# For ENGLISH: Change trait names
	if (outType == "ENGLISH") {
		names(df)[names(df) %in% c("Componentes_CHL", "Heredabilidad")] <- c("CHL_Components", "Heritability") # Rename columns "a" and "b" to "x" and "y"
		for (i in seq_len(nrow(mappings))) {
		  df[] <- lapply(df, function(col) {
			gsub(
			  mappings$SpanishName[i],
			  mappings$EnglishName[i],
			  col,
			  ignore.case = TRUE
			)
		  })
		}  
	}
 
	# Save the new table
	output_filename <- file.path (outputDir, addLabel (basename (hclTableFilename), outType))
	write.csv (df, output_filename , quote=F, row.names=F)
	return (output_filename)
}

# ---- Small helpers ----
getVarFilenameForModel <- function(model) {
  switch(model,
    "RKHS"    = "varU.dat",
    "BRR"     = "varB.dat",
    "BayesA"  = "ScaleBayesA.dat",
    "BayesB"  = "parBayesB.dat",
    "BayesC"  = "varB.dat",
    "BayesRR" = "varB.dat",
    "BL"      = "varB.dat",
    stop(sprintf("Unknown model: %s", model))
  )
}


subset_ETA <- function(ETA, keep) {
  out <- ETA
  for (nm in seq_along(ETA)) {
    if (!is.null(ETA[[nm]]$K)) out[[nm]]$K <- ETA[[nm]]$K[keep, keep, drop = FALSE]
    if (!is.null(ETA[[nm]]$X)) out[[nm]]$X <- ETA[[nm]]$X[keep, , drop = FALSE]
  }
  out
}

build_ETA <- function(geno, MODEL) {
  # Single additive component as ETA_1; kernel scaled so mean(diag(K)) = 1 for RKHS
  if (MODEL == "RKHS") {
    Z <- scale(geno, center = TRUE, scale = TRUE)
    K <- tcrossprod(Z) / ncol(Z)
    K <- K / mean(diag(K))
    return(list(list(K = K, model = "RKHS")))
  } else {
    X <- scale(geno) / sqrt(ncol(geno))
    return(list(list(X = X, model = MODEL)))
  }
}

#---------------------------------------------------------------------
# Create a table with info of NAs by trait
#---------------------------------------------------------------------
log_trait_missingness <- function(pheno, out_csv) {
  nn  <- colSums(!is.na(pheno))
  nax <- colSums(is.na(pheno))
  df  <- data.frame(trait = names(nn), n_non_missing = as.integer(nn), n_missing = as.integer(nax))
  df  <- df[order(df$n_non_missing, decreasing = TRUE),]
  write.csv(df, out_csv, row.names = FALSE)
  message("Wrote trait missingness report to: ", out_csv)
  invisible(df)
}
#--------------------------------------------------------------------
#--------------------------------------------------------------------
main ()

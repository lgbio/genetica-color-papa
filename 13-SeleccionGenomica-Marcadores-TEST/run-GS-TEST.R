#!/usr/bin/env Rscript

#------------------- For calling from a renv environment -------------------------------------------------
Sys.setenv(RENV_PROJECT = "/home/lg/BIO/agrosavia/genetica-color-papa/")
source(file.path(Sys.getenv("RENV_PROJECT"), "renv/activate.R"))
#---------------------------------------------------------------------------------------------------------

source ("lglib14.R")
library (ppGS)
library (dplyr)
library (ggplot2)
#source ("gs_plotsHCL_30.R")

#RUNTYPE = "TEST"
RUNTYPE = "FULL"

inputDir        = "inputs"
outputDir       = sprintf ("outputs-%s", RUNTYPE)
outputData      = sprintf ("data-%s", RUNTYPE)
configFile      = sprintf ("config-multi-%s.yml", RUNTYPE)
predictionsFile = sprintf ("%s/%s", outputDir, "out-CrossValidation-kfolds-GEBVs.csv")

set.seed (123)
#if (RUNTYPE == "FULL") {
#    gs_split ("inputs/genotipo.csv", "inputs/fenotipos.csv", outputDir=outputData) 
#}else {
#    gs_split ("inputs/genotipo.csv", "inputs/fenotipos.csv", outputDir=outputData, nGenos=100, nPheno=2) # For debug
#}

#--------------------------------------------------------------------
#message ("+++ Running GS for multi phenotypes...")
#predictionsKFoldsTableFilename = gs_multi (configFile, outputDir)
#setwd ("..")
#--------------------------------------------------------------------

# Parameters for columns and titles
#COLUMNS = list (prefix="Prefix", trait="HCL_component", models="Models", prediction="Predictive_ability", component="HCL_component")
#TITLES  = list (X="HCL_Components", Y="Predictive Ability", MAIN="Genomic Selection for 21 HCL Color Traits")

message ("+++ Creating plots for predictions in filename: ", predictionsFile)
gs_plotsHCL (predictionsFile, outputDir, "Predictive_ability", "GS for HCL color traits")


#!/usr/bin/Rscript
## code to prepare `DATASET` dataset goes here

# sources
sourceGeno   = read.csv ("sources/sources-genotypes-CCC-Andigena-ClusterCall2020.csv", check.names=F)
sourcePhenos = read.csv ("sources/sources-phenotypes-CCC-Andigena-ColorTraitsHCL.csv", check.names=F)

usethis::use_data(sourceGeno, overwrite=T)
usethis::use_data(sourcePhenos, overwrite=T)

# small
smallTrGn   = read.csv ("small/small-training-geno.csv", check.names=F)
smallTgGn   = read.csv ("small/small-target-geno.csv", check.names=F)
smallTrPh   = read.csv ("small/small-training-pheno.csv", check.names=F)
smallTgPh   = read.csv ("small/small-target-pheno.csv", check.names=F)
smallConfig = read.delim ("small/small-config.yml", sep="\n", header=F)

usethis::use_data(smallTrGn, overwrite=T)
usethis::use_data(smallTgGn, overwrite=T)
usethis::use_data(smallTrPh, overwrite=T)
usethis::use_data(smallTgPh, overwrite=T)
usethis::use_data(smallConfig, overwrite=T)

# full
fullTrGn    = read.csv ("full/full-training-geno.csv", check.names=F)
fullTgGn    = read.csv ("full/full-target-geno.csv", check.names=F)
fullTrPh    = read.csv ("full/full-training-pheno.csv", check.names=F)
fullTgPh    = read.csv ("full/full-target-pheno.csv", check.names=F)
fullConfig = read.delim ("full/full-config.yml", sep="\n", header=F)

usethis::use_data(fullTrGn, overwrite=T)
usethis::use_data(fullTgGn, overwrite=T)
usethis::use_data(fullTrPh, overwrite=T)
usethis::use_data(fullTgPh, overwrite=T)
usethis::use_data(fullConfig, overwrite=T)

# single
singleTrGn    = read.csv ("single/single-training-geno.csv", check.names=F)
singleTgGn    = read.csv ("single/single-target-geno.csv", check.names=F)
singleTrPh    = read.csv ("single/single-training-pheno.csv", check.names=F)
singleTgPh    = read.csv ("single/single-target-pheno.csv", check.names=F)
singleConfig  = read.delim ("single/single-config.yml", sep="\n", header=F)

usethis::use_data(singleTrGn, overwrite=T)
usethis::use_data(singleTgGn, overwrite=T)
usethis::use_data(singleTrPh, overwrite=T)
usethis::use_data(singleTgPh, overwrite=T)
usethis::use_data(singleConfig, overwrite=T)

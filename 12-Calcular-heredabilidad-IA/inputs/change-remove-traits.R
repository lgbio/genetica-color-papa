#!/usr/bin/env Rscript
library(dplyr)

# Sample data frame
source ("lglib14.R")

args = commandArgs(trailingOnly = TRUE)
inputFile = args [1]
outputFile = addLabel (inputFile, "NEWNAMES")

# Columns to remove
cols_to_remove <- c("CBaya.L", "CBaya.C", "CBaya.H", "CPulpa.L", "CPulpa.C", "CPulpa.H")
mappings  = read.csv ("traitnames-BASE.csv")

df = read.csv (inputFile)

# Remove specified columns
df <- df[ , !(names(df) %in% cols_to_remove)]

mappings$SpanishName <- as.character(mappings$SpanishName)
mappings$EnglishName <- as.character(mappings$EnglishName)
print ("...")
print ("...")
print ("...")
print ("...")


# Rename using pattern replacement
df <- df %>%
  rename_with(function(colname) {
    for (i in seq_len(nrow(mappings))) {
      colname <- gsub(mappings$SpanishName[i], mappings$EnglishName[i], colname)
    }
    colname
  })


# Keep only mapped English names present in df
cols_to_keep <- c("Registro", sort (names (df [,-1])))
print (cols_to_keep)

# Reorder df columns by that order
df <- df[, cols_to_keep, drop = FALSE]

other_cols <- setdiff(names(df), cols_to_keep)
df <- df[, c(cols_to_keep, other_cols), drop = FALSE]

view (df)


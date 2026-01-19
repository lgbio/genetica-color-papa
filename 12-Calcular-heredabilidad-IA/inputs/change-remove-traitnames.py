#!/usr/bin/env python3

import os

args = commandArgs(trailingOnly = TRUE)
inputFilename  = args [1]
outputFilename = paste0 ("NEW-", inputFilename)

traitNames = open ('traitnames-BASE.csv').readlines ()

for line in traitNames:
	values = line.split (",")
	spanish, english, description = values [0].strip(), values [1].strip(), values [2].strip()
	print (f"{spanish}, {english}, {description}")

	cmm = f'sed -i "s/{spanish}/{english}/g" {inputFilename} '
	print (cmm, "\n")

	os.system (cmm)

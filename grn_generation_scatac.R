#!/usr/bin/env Rscript

#################################################################
# scATAC preparation and GRN generation
#################################################################

library(ArchR)
library(parallel)

source("./scripts/helper_functions.R")

ArchR::addArchRThreads(threads = 16) 
ArchR::addArchRGenome("hg38")

#################################################################
# SETUP
#################################################################

scatac_folder <- "/mnt/d/scatac_input/greenleaf_23"
scatac_files <- list.files(scatac_folder, full.names = TRUE) %>%  .[grepl("\\.tsv\\.gz$", .)]
scatac_samples <- c("C_PB1", "C_PB2", "C_PB3", "C_SD1", "C_SD2", "C_SD3")

#################################################################
# CREATE ARROW FILES
#################################################################

ArrowFiles <- createArrowFiles(
  inputFiles = scatac_files,
  sampleNames = scatac_samples,
  minTSS = 4,
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

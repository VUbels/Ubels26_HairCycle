#!/usr/bin/env Rscript

#################################################################
# Description: Preproccesing data using CellBender and scDblFinder
#################################################################

####################
# 1. LIBRARY LOADING
####################

library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)
library(reticulate)

###########################
# 2. SETUP PROJECT PARAMETERS
############################

project <- "ubels26_haircycle"
input_folder <- "/mnt/d/scrna_datasets/ubels26_scrna_dataset"
preprocess_output_folder <- paste0("./preprocessing/")
scrna_data <- list.files(input_folder, recursive = FALSE, include.dirs = FALSE, pattern = ".h5")
variables <- c("anagen", "catagen", "telogen")

#########################
# 3. SETUP PY ENVIRONMENT
#########################

# Please note that for GPU support you need to manually change
# parameters in the setup_py_env() function found in helper_function.R
# Due to this being highly user dependent, questions regarding setting up
# appropriate pytorch compatibility will not be supported.
# CellBender can run without GPU support but this will take a very long time.

source("./helper_functions.R")
source("./setup_py_env.R")

py_location <- "/home/uvictor/miniconda3/bin/conda"
conda_info_env <- setup_py_env(project, py_location)
cellbender <- reticulate::import("cellbender")

#########################################
# 2. RUN CELLBENDER - AMBIENT RNA REMOVAL
#########################################

for (dataset in scrna_data) {
  
  # Construct paths
  cellbender_input_path <- file.path(input_folder, dataset)
  dataset_base <- sub("\\.h5$", "", dataset)
  
  # Create output directory
  cellbender_output_dir <- file.path(preprocess_output_folder, paste0(dataset_base, "_cellbender_results"))
  dir.create(cellbender_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Output file inside that directory
  output_name <- paste0(dataset_base, "_postcellbender.h5")
  cellbender_output_path <- file.path(cellbender_output_dir, output_name)
  
  # Build CellBender command
  cmd <- sprintf(
    "%s remove-background --input %s --output %s --cuda",
    conda_info_env$cellbender_bin,
    shQuote(cellbender_input_path),
    shQuote(cellbender_output_path)
  )

  # Run CellBender
  message("Processing: ", dataset)
  exit_code <- system(cmd)
  
  if (exit_code != 0) {
    warning("CellBender failed for: ", dataset)
  } else {
    message("Completed: ", output_name)
  }
}

################################
# 3. DATA LOADING AND INITIAL QC
################################

object.list <- list()

for (i in seq_along(objects)) {
  object <- objects[[i]]  
  stage <- variables[[i]]
  
  # Load CellBender-corrected data
  data.arna_corrected <- Read10X_h5(filename = paste0(input_folder, object), use.names = TRUE)  
  obj <- CreateSeuratObject(counts = data.arna_corrected, project = stage)
  obj$orig.ident <- stage
  
  # Calculate QC metrics
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  # QC visualization
  print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  
  plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(plot1 + plot2)
  
  object.list[[i]] <- obj
  
  rm(data.arna_corrected)
  rm(obj)
}

# 4. INITIAL QC FILTERING
for (i in seq_along(object.list)) {
  obj <- object.list[[i]]
  
  # Standard QC thresholds - adjust based on your data
  obj <- subset(obj, subset = nFeature_RNA > 100 & nFeature_RNA < 7000 & percent.mt < 30)
  
  object.list[[i]] <- obj
  cat("Remaining cells after initial QC for", unique(obj$orig.ident), "is", ncol(obj), "cells\n")
}

# 5. scDblFinder DOUBLET DETECTION
for (i in seq_along(object.list)) {
  obj <- object.list[[i]]
  stage <- unique(obj$orig.ident)
  
  cat("\n*******************************\n")
  cat("Running scDblFinder on", stage, "...\n")
  cat("********************************\n")
  
  # Convert Seurat object to SingleCellExperiment
  sce <- as.SingleCellExperiment(obj)
  
  # Run scDblFinder
  # clusters=FALSE uses the random approach (recommended for complex datasets)
  # Set seed for reproducibility
  set.seed(123)
  sce <- scDblFinder(sce, clusters = FALSE)
  
  # Extract doublet information and add to Seurat object
  obj$scDblFinder.class <- sce$scDblFinder.class
  obj$scDblFinder.score <- sce$scDblFinder.score
  
  # Print doublet statistics
  doublet_table <- table(sce$scDblFinder.class)
  cat("Doublets detected:", doublet_table["doublet"], 
      "(", round(doublet_table["doublet"]/ncol(obj)*100, 2), "%)\n")
  cat("Singlets:", doublet_table["singlet"], 
      "(", round(doublet_table["singlet"]/ncol(obj)*100, 2), "%)\n")
  
  # Update object in list
  object.list[[i]] <- obj
  
  rm(sce)
}

# 6. VISUALIZATION OF DOUBLET DETECTION
dir.create(paste0(output_folder, "scDblFinder/"), showWarnings = FALSE)

for (i in seq_along(object.list)) {
  obj <- object.list[[i]]
  stage <- unique(obj$orig.ident)
  
  # Preprocessing for visualization
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
  
  # Plot doublets on UMAP
  p1 <- DimPlot(obj, group.by = "scDblFinder.class", pt.size = 0.5) + 
    ggtitle(paste(stage, "- Doublet Classification"))
  
  p2 <- FeaturePlot(obj, features = "scDblFinder.score", pt.size = 0.5) + 
    ggtitle(paste(stage, "- Doublet Score"))
  
  # Save plots
  combined_plot <- p1 + p2
  ggsave(filename = paste0(output_folder, "scDblFinder/", stage, "_doublet_detection.png"),
         plot = combined_plot, width = 12, height = 5)
  
  print(combined_plot)
  
  object.list[[i]] <- obj
}

# 7. FILTER DOUBLETS
# Filter singlets only
object.list_filtered <- lapply(object.list, function(obj) {
  stage <- unique(obj$orig.ident)
  obj_filtered <- subset(obj, subset = scDblFinder.class == "singlet")
  cat("Remaining cells after doublet removal for", stage, ":", ncol(obj_filtered), "cells\n")
  return(obj_filtered)
})
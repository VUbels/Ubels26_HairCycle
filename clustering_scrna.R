#!/usr/bin/env Rscript

#################################################################
# Description: Clustering and evaluation of subclusters using CDI
#################################################################

# Following script contains the full clusteringpipeline for the Ubels26_HairCycle
# publication. Renv/Conda environments are dynamically set. Pytorch compatability
# has to be set by user and will not be supported.

#################################################################
# LIBRARY LOADING
#################################################################

library(Seurat)
library(SeuratDisk)

library(Nebulosa)
library(CDI)
library(reticulate)

#################################################################
# SETUP PROJECT PARAMETERS
#################################################################

project <- "ubels26_haircycle"
main_folder <- "./"
obj <- readRDS(paste0(main_folder, "post_filter_integrated_objects.RDS"))

gene_list = c("PDGFRA", "FGF7", "VIM", "EDN3", "WNT5A", "KRT1", "KRT10", "KRT6A",
              "KRT17", "KRT75", "KRT35", "KRT85", "PECAM1", "VWF", "TAGLN", "DSP",
              "MLANA", "PMEL", "SOX7","SOX9", "SOX18","SERPINA3", "PDZRN3", "FGL2",
              "CXCL14", "LGR5", "COMP", "CD34", "AQP3", "KRT79", "KRT19", "WNT3",
              "WNT10A", "WNT10B", "RGS5", "CD3D", "CD69", "TWIST2", "SOX10", "CDH9",
              "ACTA2", "EGF", "LGR6", "FGF14", "FMN2", "CTNNA2", "GPR183", "CD53",
              "CD2", "HLA-DQA2", "S100A2", "NES", "RGS5")

#################################################################
# SETUP PY ENVIRONMENT
#################################################################

# Please note that for GPU support you need to manually change
# parameters in setup_py_env.R Due to this being highly user 
# dependent, questions regarding setting up appropriate pytorch
# compatibility will not be supported. CellBender can run 
# without GPU support but this will take a very long time.

source("./scripts/helper_functions.R")
source("./scripts/setup_py_env.R")
source("./scripts/ambient_rna_removal.R")
source("./scripts/doublet_removal.R")

py_location <- "/home/uvictor/miniconda3/bin/conda"
conda_info_env <- setup_py_env(project, py_location)

#################################################################
# RUNNING BROAD MARKER GENES FOR INITIAL CLUSTERIZATION
#################################################################

# Easy visualization through Nebulosa to get a better overview of gene expression particularly
# for when low cell count has high gene expression in a particular cluster

broad_markers <- FindAllMarkers(obj, min.pct = 0.1, logfc.threshold = 0.3)

plot_marker_genes(obj = obj, 
                              genes = gene_list, 
                              cluster_col = "seurat_clusters",
                              reduction = "umap", 
                              output_dir = "./marker_genes/broad_markers", 
                              pt_size = 1,
                              outline_size = 0.25,
                              concavity = 5,
                              show_labels = TRUE,
                              eps = 2,
                              min_pts = 25,
                              outlier_percentile = 0.98)

#################################################################
# ASSIGN BROAD MARKER IDENTIFICATION TO CLUSTERS
#################################################################

broad_cluster_identification <- list(
 `0` = "Temporal.Follicle", #PDZRN3/SERPINA3/SOX9/KRT75
 `1` = "Temporal.Follicle", #PDZRN3/SOX9
 `2` = "Central.Follicle", #WNT5A/KRT19
 `3` = "Central.Follicle", #WNT5A/
 `4` = "Matrix", #KRT35/KRT85/WNT10B
 `5` = "Permanent.Follicle", #KRT1/KRT10/KRT79/AQP3
 `6` = "Permanent.Follicle", #LGR6/SOX7/WNT10A
 `7` = "Permanent.Follicle", #KRT1/KRT10/KRT79/WNT3
 `8` = "Endothelial", #PLVAP/PECAM1/VWF/SOX18/VIM
 `9` = "Pericytes", #RGS5 but also EDN3/VIM/CD34
 `10` = "Immune", #CD3D/CD53/CD69
 `11` = "Dermal.Papilla", #FGF7/TWIST2/PDGFRA
 `12` = "Melanocytes", #PMEL/MLANA/SOX10
 `13` = "Dermal.Sheath", #TALGN/VIM/ACTA2
 `14` = "Central.Follicle", #EGF/FGF14
 `15` = "Temporal.Follicle", #Likely differentiating cells? #KR17/
 `16` = "Neural.Progenitors" #CTNNA2/CDH9
)
obj$broad_cluster <- unname(unlist(broad_cluster_identification[as.character(obj$seurat_clusters)]))

visualize_percentage_clusters(seurat_obj = obj, clusters = "broad_cluster", phases = "orig.ident", output_dir = paste0(main_folder, "marker_genes"))

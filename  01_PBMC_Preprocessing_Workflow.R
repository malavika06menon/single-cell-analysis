# -------------------------------------------------------------------------
# Project: Single-Cell Analysis of PBMC 3k 
# Author: Malavika
# -------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(patchwork)

# --- Step 1: Load Data ---
# Select the 'hg19' folder from the popup window
data_path <- choose.dir()

pbmc.data <- Read10X(data.dir = data_path)

# --- Step 2: Initialize Seurat Object ---
# Min.cells = 3  -> Ignore genes that are rarely found
# Min.features = 200 -> Ignore empty droplets (dead cells)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# Check what we made
print(pbmc)

# --- Step 3: Quality Control (QC) ---
# IMPORTANT: We name this 'percent.mt' to match standard conventions
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize the quality
# Look for cells with high MT% (dead) or low nFeature (empty)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Scatter plots to see correlations
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
# Filter the data
# Rule 1: nFeature_RNA > 200  (Must have at least 200 genes detected)
# Rule 2: nFeature_RNA < 2500 (Must have fewer than 2500 genes - avoids doublets)
# Rule 3: percent.mt < 5      (Must have less than 5% mitochondrial reads - avoids dead cells)

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Check the new size of your data
pbmc
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# --- Step 5: Scaling the Data ---
# Shifts expression of each gene so that the mean is 0 and variance is 1
# This ensures highly expressed genes don't dominate the analysis
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# --- Step 6: Perform Linear Dimensionality Reduction (PCA) ---
# Compresses the data into "Principal Components" to find major sources of variation
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Visualize the PCA results
# This plot shows if your cells are starting to separate into groups
DimPlot(pbmc, reduction = "pca")

# --- Step 7: Inspect the PCA (Determine Dimensionality) ---
# We need to decide how many Principal Components (PCs) are "real" and how many are noise.

# 1. Print the top genes for the first 5 PCs
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# 2. Visualize the gene rankings (Top genes driving the separation)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

# 3. Heatmaps (The Visual Check)
# Look for clear blocks of color. If it looks like "static" TV noise, that PC is useless.
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# 4. The Elbow Plot (The Decision Maker)
# This plots the standard deviation of each PC.
# Look for the "Elbow" where the curve flattens out. That is your cutoff point.
ElbowPlot(pbmc)

# --- Step 8: Cluster the Cells ---
# We use dims = 1:10 because the Elbow Plot showed the signal fading after PC 10.
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# --- Step 9: Run UMAP (The Final Map) ---
# This creates the visual "islands" for each cell type
pbmc <- RunUMAP(pbmc, dims = 1:10)

# --- Step 10: Visualize Clusters ---
DimPlot(pbmc, reduction = "umap", label = TRUE)

# --- Step 11: Save the Processed Data ---
# Save the Seurat object so we can load it for Part 2 (Cell Type Annotation)
saveRDS(pbmc, file = "pbmc_tutorial.rds")
 



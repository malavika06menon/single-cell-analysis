 # Single-Cell RNA-Seq Analysis of PBMC 3k

## 📌 Project Overview
This project performs an end-to-end single-cell RNA sequencing analysis of 2,700 Peripheral Blood Mononuclear Cells (PBMCs) to identify distinct immune cell populations. Raw 10X Genomics data is processed, filtered, and analyzed using the **Seurat** R package to uncover transcriptionally distinct cell clusters.

## 🛠️ Tech Stack
- **Language:** R (v4.x)
- **Libraries:** Seurat, tidyverse, patchwork
- **Data Source:** 10X Genomics (PBMC 3k dataset)

## 📂 Repository Structure
- **`01_PBMC_Preprocessing_Workflow.R`**  
  Core analysis script implementing the full preprocessing and initial analysis pipeline:
  - **Quality Control:** Removal of low-quality cells based on feature counts and mitochondrial RNA content (<5%)
  - **Normalization:** Global-scaling normalization using LogNormalize
  - **Feature Selection:** Identification of highly variable genes
  - **Dimensionality Reduction:** Principal Component Analysis (PCA)
  - **Clustering:** Graph-based clustering
  - **Visualization:** UMAP projection for cluster visualization

## 📊 Key Results
The analysis identified **9 distinct transcriptional clusters**, corresponding to major immune cell types, including:
- **T cells** (IL7R⁺)
- **B cells** (MS4A1⁺)
- **Monocytes** (CD14⁺, FCGR3A⁺)
- **Natural Killer (NK) cells** (GNLY⁺)
- **Platelets** (PPBP⁺)

Cluster identities were inferred based on established canonical marker genes.

## 🚀 How to Run
1. Clone this repository.
2. Download the PBMC 3k dataset from 10X Genomics.
3. Open and run `01_PBMC_Preprocessing_Workflow.R` in RStudio.

## 📈 Next Steps
- Differential expression analysis to identify cluster-specific marker genes
- Cell type annotation validation
- Extension to disease vs. healthy datasets for biological interpretation

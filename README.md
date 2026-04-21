# Singlecell-analysis-using-publicly-available-datasets
# 🫁 Single-cell RNA-seq Analysis of Stage III Squamous Cell Lung Carcinoma (NSCLC)

This repository contains analysis of a **10x Genomics single-cell RNA-seq dataset** derived from Stage III squamous cell lung carcinoma (NSCLC) tumor samples.
The dataset consists of dissociated tumor cells (DTCs) processed using the **10x Genomics Chromium X platform (5' Gene Expression v2 chemistry)**.

---

# 🧬 Dataset Overview

## 🫁 Biological source
- Disease: Non-small cell lung cancer (NSCLC)
- Subtype: Squamous cell carcinoma
- Stage: III
- Sample type: Dissociated tumor cells (DTCs)
- Species: Human
- Tissue: Lung
- Preservation: Cryopreserved

---

## 🔬 Single-cell RNA-seq data

- Platform: 10x Genomics Chromium X
- Chemistry: 5' Gene Expression (v2)
- Input cells: ~4,000
- Recovered cells: 2,616
- Sequencing depth: ~25,000 reads per cell
- Library type: Paired-end Illumina sequencing

### Read structure:
- Read 1: Cell barcode + UMI
- i5 / i7: Sample indexes
- Read 2: Transcript sequence

---

# 📦 Data files used
Only the **10x Genomics gene expression output** is used:
- `filtered_feature_bc_matrix.h5` (recommended)
  OR
- `filtered_feature_bc_matrix/` (MTX format)

---

# 🧠 Objectives

This project focuses on:
- 🧬 Quality control of single-cell RNA-seq data
- 🧬 Normalization and scaling
- 🧬 Dimensionality reduction (PCA, UMAP)
- 🧬 Cell clustering
- 🧬 Identification of cell types in NSCLC tumor
- 🧬 Marker gene analysis (tumor vs immune cells)

---

# ⚙️ Analysis workflow

## 1. Data loading
- Import 10x Genomics filtered matrix

## 2. Quality control
- Filter cells based on:
  - gene count
  - UMI count
  - mitochondrial percentage

## 3. Normalization
- Log normalization or SCTransform

## 4. Feature selection
- Highly variable genes (HVGs)

## 5. Dimensionality reduction
- PCA
- UMAP / t-SNE

## 6. Clustering
- Leiden / Louvain clustering

## 7. Cell annotation
- Identify:
  - tumor cells
  - immune cells (T, B, NK, macrophages)
  - stromal cells

## 8. Differential expression
- Marker gene analysis per cluster

---

# 🧰 Tools used

### R workflow
- Seurat

### Python workflow
- Scanpy

---

# 🚫 What is NOT included

- No Oxford Nanopore long-read analysis
- No isoform / splicing analysis
- No multi-omics integration
You may explore the additional data from https://www.10xgenomics.com/datasets/3k-human-squamous-cell-lung-carcinoma-dtcs-chromium-x-2-standard

---

# 📊 Expected outputs

- UMAP plots of cell clusters
- Heatmaps of marker genes
- Cluster annotations
- Tumor vs immune composition analysis

---

# 📜 Metadata

| Attribute | Value |
|----------|------|
| Species | Homo sapiens |
| Tissue | Lung |
| Disease | NSCLC (Squamous cell carcinoma) |
| Platform | 10x Genomics Chromium X |
| Chemistry | 5' Gene Expression v2 |
| Cells recovered | 2,616 |

---

# 🚀 Purpose

This project provides a **standard single-cell RNA-seq analysis pipeline for lung cancer tumor cells**, focusing on tumor heterogeneity and microenvironment composition.

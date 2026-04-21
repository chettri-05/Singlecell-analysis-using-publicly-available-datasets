
## 🫁 Single-cell RNA-seq Analysis of NSCLC Tumor Cells (Demonstration Project)

This repository contains a **demonstration workflow for single-cell RNA sequencing (scRNA-seq) analysis** using a publicly available dataset of Stage III squamous cell lung carcinoma (NSCLC) tumor cells.

> ⚠️ **Disclaimer:**  
> This project is intended **for educational and demonstration purposes only**. It is not a clinical study and does not claim any novel biological discovery.

---

## 🧬 Dataset Overview

### 🫁 Biological context
- Disease: Non-small cell lung cancer (NSCLC)
- Subtype: Squamous cell carcinoma
- Stage: III
- Sample type: Dissociated tumor cells (DTCs)
- Species: Human
- Tissue: Lung
- Preservation: Cryopreserved

---

### 🔬 Single-cell RNA-seq data

- Platform: 10x Genomics Chromium X
- Chemistry: 5' Gene Expression (v2)
- Input cells: ~4,000
- Recovered cells: 2,616
- Sequencing depth: ~25,000 reads per cell

---

## 🎯 Project Objective & Structure

This repository demonstrates a **standard scRNA-seq analysis pipeline**, including:

- Quality control (QC) of single-cell data
- Normalization and scaling
- Identification of highly variable genes
- Dimensionality reduction (PCA, UMAP)
- Unsupervised clustering
- Cell type annotation
- Marker gene identification

```
NSCLC-scRNA-demo/
│
├── README.md
├── LICENSE
├── .gitignore
│
├── data/
│   ├── raw/                  # downloaded 10x files (optional)
│   └── processed/            # filtered_feature_bc_matrix.h5
│
├── scripts/
│   ├── 01_load_qc.R
│   ├── 02_normalization.R
│   ├── 03_dimensionality_reduction.R
│   ├── 04_clustering.R
│   ├── 05_marker_analysis.R
│
│   ├── 01_load_qc.py
│   ├── 02_preprocessing.py
│   ├── 03_umap_clustering.py
│
├── results/
│   ├── figures/
│   ├── tables/
│
├── notebooks/
│   ├── seurat_analysis.ipynb
│   ├── scanpy_analysis.ipynb
│
└── docs/
    └── workflow_diagram.png
```
---

## ⚙️ Analysis Workflow

### 1. Data loading
- Import 10x Genomics filtered gene expression matrix

### 2. Quality control
- Filter low-quality cells based on:
  - gene count thresholds
  - UMI counts
  - mitochondrial gene percentage

### 3. Data normalization
- Log normalization or SCTransform

### 4. Feature selection
- Identification of highly variable genes (HVGs)

### 5. Dimensionality reduction
- Principal Component Analysis (PCA)
- UMAP visualization

### 6. Clustering
- Leiden / Louvain clustering algorithm

### 7. Cell type annotation
- Identification of major cell populations:
  - tumor epithelial cells
  - immune cells (T cells, B cells, NK cells, macrophages)
  - stromal cells

### 8. Differential expression analysis
- Identification of cluster-specific marker genes

---

## 🧰 Tools & Technologies

### R-based analysis
- Seurat

### Python-based analysis
- Scanpy
- AnnData framework

### Visualization
- UMAP / t-SNE plots
- Heatmaps
- Violin plots

---

## 📦 Input Data

Only processed 10x Genomics gene expression data is used:

- `filtered_feature_bc_matrix.h5` (preferred)
OR
- `filtered_feature_bc_matrix/` directory

---

## 📊 Expected Outputs

- Cell clustering visualization (UMAP)
- Cluster annotation maps
- Marker gene expression heatmaps
- Tumor vs immune composition analysis

---

## 🚫 Scope Limitations

This repository does NOT include:
- Long-read sequencing analysis (Oxford Nanopore)
- Isoform or splicing-level analysis
- One can add further analysis taking the data from https://www.10xgenomics.com/datasets/3k-human-squamous-cell-lung-carcinoma-dtcs-chromium-x-2-standard

---

## 🧠 Purpose

This project is designed as a:

> 🧬 **Educational demonstration of standard single-cell RNA-seq analysis workflows in lung cancer research**

It is intended for:
- students
- beginners in bioinformatics
- scRNA-seq pipeline demonstration
- methodological learning

---

## 📜 Data Source

Publicly available 10x Genomics dataset of NSCLC tumor cells.  
All data used complies with the original data usage and licensing terms.

---

## ⚖️ License

This repository is for educational purposes only.  
No clinical or diagnostic claims are made.

---

# 🚀 Author Note

This project demonstrates a reproducible workflow for scRNA-seq analysis using standard bioinformatics tools and best practices.

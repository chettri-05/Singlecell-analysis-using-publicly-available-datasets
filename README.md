
## 🫁 scRNA-seq Analysis: 3k Human Squamous Cell Lung Carcinoma Disseminated Tumor Cells (DTCs) (Demonstration Project)

![R](https://img.shields.io/badge/R-%3E%3D4.3-blue) ![Seurat](https://img.shields.io/badge/Seurat-v5-green) ![License](https://img.shields.io/badge/license-MIT-lightgrey)

End-to-end single-cell RNA-seq analysis pipeline for the **3k Human Squamous Cell Lung Carcinoma DTC** dataset (10x Chromium 5' v2), built with [Seurat v5](https://satijalab.org/seurat/) and standard Bioconductor tools.
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

## Repository Structure

```
.
├── analysis.R                        # Full analysis script (Sections 1–14)
├── README.md
├── output/
│   ├── lung_dtc_final.rds            # Final Seurat object
│   ├── lung_dtc_all_markers.csv      # Marker genes for every cluster
│   ├── epithelial_markers.csv        # Markers for epithelial subclusters
│   ├── cluster_celltype_composition.csv
│   └── images/
│       └── lung_dtc_umap_final.png   # Publication-quality UMAP
```

---

## ⚙️ Pipeline Overview

### 1. Download input files, Environment Setup & Package Installation
📦 Input Data: Only processed 10x Genomics gene expression data is used:
- `filtered_feature_bc_matrix.h5` (preferred)
OR
- `filtered_feature_bc_matrix/` directory

📦 Install all required CRAN, Bioconductor, and GitHub packages. Includes Seurat v5, SCTransform, SingleR, AUCell, ComplexHeatmap, and msigdbr.

### 2. Load Libraries
Load all libraries needed for the full pipeline.

### 3. Load Data & Create Seurat Object
Read the filtered HDF5 file with `Read10X_h5()` and initialise a Seurat object. Genes expressed in fewer than 3 cells and cells with fewer than 200 detected genes are discarded at this stage.

### 4. Quality Control (QC)
- Compute **mitochondrial %** (`percent.mt`) and **ribosomal %** (`percent.rb`)
- Visualise distributions with violin plots and scatter plots
- Filter cells: `nFeature_RNA` 200–6,000, `percent.mt` < 20%
- Re-visualise after filtering to confirm thresholds are appropriate

### 5. Normalization with SCTransform
Replace the classical `NormalizeData → FindVariableFeatures → ScaleData` workflow with a single `SCTransform()` call. Mitochondrial % is regressed out automatically. Identifies 3,000 variable features by default.

### 6. Dimensionality Reduction: PCA
- Run PCA on SCTransform variable features
- Assess dimensionality with `ElbowPlot()` (variance explained per PC)
- Inspect top loading genes per PC and dimensional heatmaps

### 7. Clustering
- Build a KNN graph with `FindNeighbors()` (dims 1:20)
- Cluster with the Louvain algorithm via `FindClusters()` (resolution = 0.5)
- Inspect cluster sizes with `table(Idents())`

### 8. Non-linear Dimensionality Reduction: UMAP & t-SNE
- Run `RunUMAP()` and `RunTSNE()` (dims 1:20)
- Compare UMAP and t-SNE layouts side-by-side

### 9. Marker Gene Identification
- `FindAllMarkers()` (SCT assay, only positive markers, `logfc.threshold = 0.25`)
- Top 10 marker genes per cluster exported as CSV
- Known lung cancer and immune markers visualised as `FeaturePlot` and `DotPlot`
- Top 5 per cluster summarised in a `DoHeatmap`

### 10. Automated Cell Type Annotation with SingleR
- Reference datasets: **Blueprint/ENCODE** (immune focused) and **Human Primary Cell Atlas** (broad coverage)
- Consensus label: Blueprint label where available, HPCA otherwise
- Annotation confidence inspected with `plotScoreHeatmap()`
- UMAP coloured by each annotation and the consensus

### 11. AUCell: Gene Set Activity Scoring
- Custom gene sets covering 10 biologically relevant programs (Epithelial Tumor, EMT, T cells, CD8 Cytotoxic, CD4 Helper, NK, B cells, Myeloid, Fibroblast, Endothelial)
- Rankings built from SCT log-normalised matrix
- AUC scores stored in metadata; dominant program assigned per cell
- Bimodal threshold exploration with `AUCell_exploreThresholds()`

### 12. MSigDB Hallmark Pathway Scoring
- All 50 Hallmark pathways downloaded via `msigdbr`
- AUC scores calculated per cell per pathway
- Key cancer pathways visualised: **IFN-γ Response**, **Hypoxia**, **EMT**, **MYC Targets**, **TNFα–NFκB**, **G2M Checkpoint**
- Violin plots per cell type for selected pathways

### 13. Subset Analysis — Tumor Epithelial & EMT Cells
- Subset cells annotated as epithelial/keratinocyte
- Re-run SCTransform, PCA, clustering, and UMAP on the subset
- Find within-subtype marker genes
- Export markers and heatmap

### 14. Save Results
- Final Seurat object saved as `.rds`
- Marker tables, composition tables, and UMAP image saved to `output/`

---

## Requirements

**R ≥ 4.3.0**

| Package | Source | Purpose |
|---|---|---|
| Seurat ≥ 5.0 | CRAN | Core single-cell framework |
| SeuratObject | CRAN | Seurat S4 object |
| sctransform | CRAN | SCTransform normalisation |
| BPCells | r-universe | High-performance matrix ops |
| presto | r-universe | Fast Wilcoxon DE test |
| glmGamPoi | r-universe | GLM for SCTransform |
| SingleR | Bioconductor | Automated annotation |
| celldex | Bioconductor | Reference datasets for SingleR |
| AUCell | Bioconductor | Gene set activity scoring |
| ComplexHeatmap | Bioconductor | Heatmap visualisation |
| SingleCellExperiment | Bioconductor | SCE container |
| msigdbr | CRAN | MSigDB gene sets |
| GSEABase | Bioconductor | Gene set objects |
| ggplot2 | CRAN | Plotting |
| patchwork | CRAN | Plot composition |
| dplyr | CRAN | Data wrangling |
| viridis | CRAN | Colour palettes |

---

## Quick Start

```r
# 1. Clone the repository
# git clone https://github.com/<your-username>/lung-dtc-scrna.git

# 2. Place the input HDF5 file in the project root:
#    SC5pv2_GEX_Human_Lung_Carcinoma_DTC_filtered_feature_bc_matrix.h5

# 3. Open analysis.R and run from top to bottom (or source it):
source("analysis.R")

# Outputs will appear in output/ and output/images/
```
---
```
1. Environment Setup & Package Installation
r# Install BiocManager
install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("ComplexHeatmap", "celldex", "SingleR", "AUCell"), 
                     ask = FALSE, update = FALSE)

# Install CRAN packages
install.packages(c(
  "ggplot2", "ggrepel", "RColorBrewer", "patchwork",
  "dplyr", "tibble", "viridis", "Seurat", "pheatmap",
  "msigdbr", "gridExtra"
))

# Install sctransform
install.packages("sctransform")

# Install BPCells, presto, glmGamPoi from universe repos
setRepositories(ind = 1:3, addURLs = c(
  'https://satijalab.r-universe.dev',
  'https://bnprks.r-universe.dev/'
))
install.packages(c("BPCells", "presto", "glmGamPoi"))

# Install GitHub packages
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("satijalab/seurat-data", quiet = TRUE)
remotes::install_github("satijalab/azimuth", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)
2. Load Libraries
rlibrary(Seurat)
library(SeuratObject)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(tibble)
library(BPCells)
library(presto)
library(glmGamPoi)
library(ComplexHeatmap)
library(SingleCellExperiment)
library(SingleR)
library(viridis)
library(celldex)
library(AUCell)
library(pheatmap)
library(sctransform)
library(msigdbr)
library(GSEABase)
library(gridExtra)
3. Load Data & Create Seurat Object
r# Load the HDF5 file (Chromium 5' v2 output)
data <- Read10X_h5("SC5pv2_GEX_Human_Lung_Carcinoma_DTC_filtered_feature_bc_matrix.h5")

# Create Seurat object
# 36601 features x 2616 cells
seurat_obj <- CreateSeuratObject(
  counts  = data,
  project = "LungCarcinomaDTC",
  min.cells    = 3,    # keep genes detected in ≥ 3 cells
  min.features = 200   # keep cells with ≥ 200 genes detected
)

seurat_obj
dim(seurat_obj)
4. Quality Control (QC)
r# --- Compute QC metrics ---
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

# Inspect metadata
head(seurat_obj@meta.data, 5)

# --- Visualize before filtering ---
VlnPlot(seurat_obj,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
        ncol = 4, pt.size = 0.1) +
  ggtitle("QC Metrics Before Filtering")

# Scatter plots: count vs mt% and count vs feature count
p1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  ggtitle("UMI Count vs Mitochondrial %")
p2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  ggtitle("UMI Count vs Gene Count")
p3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.rb") +
  ggtitle("UMI Count vs Ribosomal %")
p1 + p2 + p3

# --- Filter cells ---
# Thresholds chosen for a ~2600-cell tumor DTC dataset:
#   nFeature_RNA > 200       → remove empty droplets
#   nFeature_RNA < 6000      → remove likely doublets (tumor cells express more genes)
#   percent.mt   < 20        → tumor datasets tolerate slightly higher mt% than PBMCs
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > 200 &
           nFeature_RNA < 6000 &
           percent.mt < 20
)

cat("Cells after QC filtering:", ncol(seurat_obj), "\n")

# --- Visualize after filtering ---
VlnPlot(seurat_obj,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
        ncol = 4, pt.size = 0.1) +
  ggtitle("QC Metrics After Filtering")
5. Normalization with SCTransform
r# SCTransform replaces NormalizeData + FindVariableFeatures + ScaleData
# Automatically regresses out sequencing depth; optionally also regress mt%
seurat_obj <- SCTransform(
  seurat_obj,
  vars.to.regress = "percent.mt",  # regress out mitochondrial contamination
  verbose         = FALSE
)

# Inspect top variable features
top10_var <- head(VariableFeatures(seurat_obj), 10)
cat("Top 10 variable features:\n")
print(top10_var)

# Plot variable features
p_var1 <- VariableFeaturePlot(seurat_obj)
p_var2 <- LabelPoints(plot = p_var1, points = top10_var, repel = TRUE)
p_var1 + p_var2
6. Dimensionality Reduction: PCA
r# Run PCA on SCT variable features
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

# Elbow plot: choose number of PCs
ElbowPlot(seurat_obj, ndims = 40) +
  ggtitle("Elbow Plot: Variance Explained per PC")

# Inspect top genes driving each PC
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)

# PCA plots
VizDimLoadings(seurat_obj, dims = 1:4, reduction = "pca")
DimPlot(seurat_obj, reduction = "pca") + ggtitle("PCA Plot")

# Dimensional heatmaps
DimHeatmap(seurat_obj, dims = 1:15, cells = 500, balanced = TRUE)
7. Clustering
r# Build KNN graph (use dims informed by elbow plot; adjust as needed)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)

# Louvain clustering — resolution 0.5 is a good start for ~2600 cells
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

cat("Number of clusters found:", length(levels(Idents(seurat_obj))), "\n")
head(Idents(seurat_obj), 5)
table(Idents(seurat_obj))
8. Non-linear Dimensionality Reduction: UMAP & t-SNE
r# UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("UMAP — Louvain Clusters (Lung Carcinoma DTC)")

# t-SNE (alternative view)
seurat_obj <- RunTSNE(seurat_obj, dims = 1:20, verbose = FALSE)
DimPlot(seurat_obj, reduction = "tsne", label = TRUE, repel = TRUE) +
  ggtitle("t-SNE — Louvain Clusters (Lung Carcinoma DTC)")

# Side-by-side
p_umap <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle("UMAP")
p_tsne <- DimPlot(seurat_obj, reduction = "tsne", label = TRUE) + ggtitle("t-SNE")
p_umap + p_tsne
9. Marker Gene Identification
r# --- Markers for each cluster vs all others ---
seurat_obj.markers <- FindAllMarkers(
  seurat_obj,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  assay           = "SCT"
)

# Top 10 markers per cluster (by log2FC)
top10_markers <- seurat_obj.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

print(top10_markers)

# Save markers table
write.csv(seurat_obj.markers, file = "output/lung_dtc_all_markers.csv", row.names = TRUE)

# --- Heatmap of top 5 markers per cluster ---
top5_markers <- seurat_obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup()

DoHeatmap(seurat_obj, features = top5_markers$gene, assay = "SCT") +
  ggtitle("Top 5 Marker Genes per Cluster") + NoLegend()

# --- Known lung cancer / immune marker genes ---
known_markers <- c(
  "EPCAM", "KRT19", "KRT18",   # Epithelial / tumor cells
  "VIM", "CDH2", "FN1",        # Mesenchymal / EMT
  "CD3D", "CD3E",              # T cells
  "MS4A1", "CD79A",            # B cells
  "LYZ", "CD14", "FCGR3A",    # Myeloid / Monocytes
  "NKG7", "GNLY",              # NK cells
  "PECAM1", "VWF",             # Endothelial
  "COL1A1", "FAP"              # Fibroblasts
)

FeaturePlot(
  seurat_obj,
  features  = known_markers,
  reduction = "umap",
  ncol      = 4,
  min.cutoff = "q10",
  max.cutoff = "q90"
) + ggtitle("Known Lung Cancer & Immune Markers")

DotPlot(seurat_obj, features = known_markers, assay = "SCT") +
  RotatedAxis() + ggtitle("Marker Expression Across Clusters")
10. Automated Cell Type Annotation with SingleR
r# Load reference datasets
ref_blueprint <- celldex::BlueprintEncodeData()
ref_hpca      <- celldex::HumanPrimaryCellAtlasData()

# Convert Seurat to SingleCellExperiment
sce_lung <- as.SingleCellExperiment(seurat_obj)

# Run SingleR with Blueprint/ENCODE (good for immune + cancer cell types)
singleR_bp <- SingleR(
  test            = sce_lung,
  assay.type.test = "logcounts",
  ref             = ref_blueprint,
  labels          = ref_blueprint$label.main
)

# Run SingleR with HPCA
singleR_hpca <- SingleR(
  test            = sce_lung,
  assay.type.test = "logcounts",
  ref             = ref_hpca,
  labels          = ref_hpca$label.main
)

# Consensus: use Blueprint labels where available; fall back to HPCA
results_consensus <- ifelse(
  !is.na(singleR_bp$pruned.labels),
  singleR_bp$pruned.labels,
  singleR_hpca$pruned.labels
)

# Store annotations in Seurat metadata
seurat_obj@meta.data$celltype_blueprint <- singleR_bp$pruned.labels
seurat_obj@meta.data$celltype_hpca      <- singleR_hpca$pruned.labels
seurat_obj@meta.data$celltype_consensus <- results_consensus

head(seurat_obj@meta.data[, c("celltype_blueprint", "celltype_hpca", "celltype_consensus")])

# --- Annotation score heatmaps ---
plotScoreHeatmap(singleR_bp,  main = "SingleR Scores — Blueprint/ENCODE")
plotScoreHeatmap(singleR_hpca, main = "SingleR Scores — HPCA")

# --- UMAP colored by annotation ---
DimPlot(seurat_obj, reduction = "umap", group.by = "celltype_blueprint",
        label = TRUE, repel = TRUE) +
  ggtitle("Blueprint/ENCODE Annotation")

DimPlot(seurat_obj, reduction = "umap", group.by = "celltype_hpca",
        label = TRUE, repel = TRUE) +
  ggtitle("HPCA Annotation")

DimPlot(seurat_obj, reduction = "umap", group.by = "celltype_consensus",
        label = TRUE, repel = TRUE) +
  ggtitle("Consensus Cell Type Annotation — Lung Carcinoma DTC")

# Cell type distribution
table(seurat_obj@meta.data$celltype_consensus)
11. AUCell: Gene Set Activity Scoring
r# --- Define lung cancer–relevant gene sets ---
markers.lung <- list(
  Epithelial_Tumor = c("EPCAM", "KRT19", "KRT18", "KRT8", "CLDN4"),
  EMT              = c("VIM", "CDH2", "FN1", "SNAI1", "TWIST1", "ZEB1"),
  T_cells          = c("CD3D", "CD3E", "CD3G", "TRAC"),
  CD8_Cytotoxic    = c("CD8A", "CD8B", "GZMB", "PRF1", "NKG7"),
  CD4_Helper       = c("CD4", "IL7R", "CCR7", "LTB"),
  NK_cells         = c("NKG7", "GNLY", "KLRD1", "NCR1"),
  B_cells          = c("MS4A1", "CD79A", "CD74"),
  Myeloid          = c("LYZ", "CD14", "S100A8", "S100A9", "CTSS"),
  Fibroblast       = c("COL1A1", "COL3A1", "FAP", "ACTA2"),
  Endothelial      = c("PECAM1", "VWF", "CDH5", "ENG")
)

# Build GeneSetCollection
all_sets <- GeneSetCollection(
  lapply(names(markers.lung), function(x) {
    GeneSet(markers.lung[[x]], setName = x)
  })
)

# Extract log-normalised expression matrix from SCT assay
sce_auc <- as.SingleCellExperiment(seurat_obj)
exprMat  <- logcounts(sce_auc)

# Build gene rankings
rankings <- AUCell_buildRankings(exprMat, plotStats = FALSE, verbose = FALSE)

# Calculate AUC scores
cell.aucs <- AUCell_calcAUC(all_sets, rankings)

# Extract AUC matrix and assign dominant label per cell
auc.mat <- t(getAUC(cell.aucs))
seurat_obj$AUCell_label <- colnames(auc.mat)[max.col(auc.mat)]

table(seurat_obj$AUCell_label)

# UMAP colored by AUCell label
DimPlot(seurat_obj, group.by = "AUCell_label", label = TRUE, repel = TRUE) +
  ggtitle("AUCell Gene Set Activity — Lung Carcinoma DTC")

# Threshold distributions (check bimodality)
AUCell_exploreThresholds(cell.aucs, plotHist = TRUE, assign = TRUE)

# Add individual AUC scores as metadata for FeaturePlot
seurat_obj <- AddMetaData(seurat_obj, as.data.frame(auc.mat))

FeaturePlot(
  seurat_obj,
  features  = c("Epithelial_Tumor", "EMT", "T_cells", "Myeloid"),
  reduction = "umap",
  ncol      = 2,
  min.cutoff = "q10",
  max.cutoff = "q90"
) + ggtitle("AUCell Scores — Key Cell Programs")
12. MSigDB Hallmark Pathway Scoring
r# Download Hallmark gene sets (H collection)
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)

hallmark_gs <- GeneSetCollection(
  lapply(names(hallmark_list), function(x) {
    GeneSet(unique(hallmark_list[[x]]), setName = x)
  })
)

cat("Number of Hallmark pathways:", length(hallmark_gs), "\n")

# Re-use rankings from Section 11 (or rebuild if needed)
hallmark.aucs   <- AUCell_calcAUC(hallmark_gs, rankings)
hallmark.mat    <- t(assay(hallmark.aucs))

# Add Hallmark AUC scores to Seurat metadata
seurat_obj <- AddMetaData(seurat_obj, as.data.frame(hallmark.mat))

# --- Feature plots: cancer-relevant Hallmark pathways ---
p_ifn  <- FeaturePlot(seurat_obj, features = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                      cols = c("grey90", "darkred"),   min.cutoff = "q10", max.cutoff = "q90") +
           ggtitle("IFN-γ Response")

p_hyp  <- FeaturePlot(seurat_obj, features = "HALLMARK_HYPOXIA",
                      cols = c("grey90", "darkorange"), min.cutoff = "q10", max.cutoff = "q90") +
           ggtitle("Hypoxia")

p_emts <- FeaturePlot(seurat_obj, features = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                      cols = c("grey90", "steelblue"),  min.cutoff = "q10", max.cutoff = "q90") +
           ggtitle("EMT (Hallmark)")

p_myc  <- FeaturePlot(seurat_obj, features = "HALLMARK_MYC_TARGETS_V1",
                      cols = c("grey90", "purple"),     min.cutoff = "q10", max.cutoff = "q90") +
           ggtitle("MYC Targets V1")

p_tnf  <- FeaturePlot(seurat_obj, features = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                      cols = c("grey90", "forestgreen"), min.cutoff = "q10", max.cutoff = "q90") +
           ggtitle("TNFα–NFκB")

p_pro  <- FeaturePlot(seurat_obj, features = "HALLMARK_G2M_CHECKPOINT",
                      cols = c("grey90", "firebrick"),  min.cutoff = "q10", max.cutoff = "q90") +
           ggtitle("G2M Checkpoint (Proliferation)")

(p_ifn | p_hyp | p_emts) / (p_myc | p_tnf | p_pro)

# Violin plots: pathway activity per cell type
VlnPlot(
  seurat_obj,
  features  = c(
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_HYPOXIA",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE"
  ),
  group.by  = "celltype_consensus",
  pt.size   = 0,
  ncol      = 1
)
13. Subset Analysis — Tumor Epithelial & EMT Cells
r# Set identity to consensus annotation
seurat_obj <- SetIdent(seurat_obj, value = "celltype_consensus")
levels(seurat_obj)

# Subset: epithelial / carcinoma-like cells (adjust names to match your annotation output)
epithelial_ids <- grep("Epithelial|epithelial|Keratinocyte|keratinocyte",
                        levels(seurat_obj), value = TRUE)

if (length(epithelial_ids) > 0) {
  epi_sub <- subset(seurat_obj, idents = epithelial_ids)

  DefaultAssay(epi_sub) <- "RNA"
  epi_sub <- SCTransform(epi_sub, vars.to.regress = "percent.mt", verbose = FALSE)
  epi_sub <- RunPCA(epi_sub)
  epi_sub <- FindNeighbors(epi_sub, dims = 1:15)
  epi_sub <- FindClusters(epi_sub, resolution = 0.4)
  epi_sub <- RunUMAP(epi_sub, dims = 1:15)

  DimPlot(epi_sub, label = TRUE) + ggtitle("Epithelial/Tumor Subclusters")

  epi.markers <- FindAllMarkers(epi_sub, only.pos = TRUE,
                                min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")

  epi_top5 <- epi.markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)
  write.csv(epi.markers, "output/epithelial_markers.csv", row.names = FALSE)

  DoHeatmap(epi_sub, features = epi_top5$gene, assay = "SCT") +
    NoLegend() + ggtitle("Top Markers — Epithelial Subclusters")
} else {
  message("No epithelial clusters found with current annotation — adjust `epithelial_ids` grep pattern.")
}
14. Save Results
r# Create output directory if not present
dir.create("output", showWarnings = FALSE)
dir.create("output/images", showWarnings = FALSE)

# Save final Seurat object
saveRDS(seurat_obj, file = "output/lung_dtc_final.rds")

# Final UMAP — publication-quality
final_umap <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by  = "celltype_consensus",
  label     = TRUE,
  repel     = TRUE,
  label.size = 4,
  pt.size   = 0.5
) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  ggtitle("3k Human Lung Carcinoma DTCs — Cell Type Annotation") +
  theme(
    axis.title      = element_text(size = 14),
    legend.text     = element_text(size = 10),
    plot.title      = element_text(size = 16, face = "bold")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 6)))

ggsave("output/images/lung_dtc_umap_final.png",
       plot = final_umap, width = 12, height = 8, dpi = 300)

# Cluster composition table
cluster_summary <- as.data.frame(table(
  Cluster   = Idents(seurat_obj),
  CellType  = seurat_obj$celltype_consensus
))
write.csv(cluster_summary, "output/cluster_celltype_composition.csv", row.names = FALSE)

cat("Analysis complete. Files saved to output/\n")

```

## Key Results

| Step | Output |
|---|---|
| Post-QC cells | ~2,400–2,550 (dataset-dependent) |
| PCs used | 20 |
| Clusters (res 0.5) | ~8–12 |
| Cell types (SingleR consensus) | Epithelial, T cells, NK, B cells, Myeloid, Fibroblasts, Endothelial |
| Top Hallmark programs | EMT, Hypoxia, IFN-γ Response, MYC Targets |

---

## References

- [Seurat v5](https://satijalab.org/seurat/) — Hao et al. (2023) *Nature Methods*
- [SCTransform](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1) — Hafemeister & Satija (2019) *Genome Biology*
- [SingleR](https://bioconductor.org/packages/SingleR/) — Aran et al. (2019) *Nature Immunology*
- [AUCell](https://bioconductor.org/packages/AUCell/) — Aibar et al. (2017) *Nature Methods*
- [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb) — Liberzon et al. (2015) *Cell Systems*
- [10x Genomics Lung Carcinoma DTC Dataset](https://www.10xgenomics.com/datasets)

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

MIT License. See `LICENSE` for details.
---

# 🚀 Author Note

This project demonstrates a reproducible workflow for scRNA-seq analysis using standard bioinformatics tools and best practices.

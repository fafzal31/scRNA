# 🔭 Scanpy scRNA-seq Analysis
### Raw Count Matrix → Annotated Cell Type Clusters

> Complete downstream single-cell RNA-seq analysis of the **PBMC 3k dataset** using **Scanpy** in Python — from count loading through quality control, dimensionality reduction, clustering, and final cell type annotation.

---

## Goal

**Raw count matrix → Biologically meaningful, annotated cell type clusters with visualizations**

---

## 📁 Folder Contents

```
scanpy_analysis/
│
├── notebooks/
│   └── pbmc3k_analysis.py          # Full analysis pipeline
│
├── figures/
│   ├── 01_qc_violin.png            # QC metrics violin plot
│   ├── 02_qc_scatter.png           # Total counts vs % mitochondrial reads
│   ├── 03_highly_variable.png      # HVG dispersion plot
│   ├── 04_pca_variance.png         # PCA variance ratio (elbow plot)
│   ├── 05_umap_markers.png         # UMAP colored by marker gene expression
│   ├── 06_umap_leiden.png          # UMAP colored by Leiden cluster
│   ├── 07_dotplot.png              # Marker gene dotplot across clusters
│   └── 08_umap_annotated.png       # Final annotated UMAP
│
└── README.md
```

---

## 📊 Dataset

| Property | Details |
|---|---|
| **Sample** | 3k PBMCs from a Healthy Donor |
| **Source** | `sc.datasets.pbmc3k()` — 10x Genomics |
| **Cells (raw)** | 2,700 |
| **Genes (raw)** | 32,738 |
| **Cells after QC** | **2,638** |
| **HVGs selected** | **~1,838** |

---

## Analysis Pipeline

```
sc.datasets.pbmc3k()
        │
        ▼  Filter cells (<200 genes) and genes (<3 cells)
        │
        ▼  Flag mitochondrial genes → calculate QC metrics
        │
        ▼  Remove low-quality cells (>2500 genes or >5% MT)
        │
        ▼  Normalize to 10,000 counts/cell → log1p transform
        │
        ▼  Select highly variable genes (~1,838)
        │
        ▼  Regress out total_counts + pct_counts_mt
        │
        ▼  Scale to zero mean, unit variance (max_value=10)
        │
        ▼  PCA — 40 components (svd_solver='arpack')
        │
        ▼  KNN neighbor graph (k=10, n_pcs=40)
        │
        ▼  Leiden clustering (resolution=0.9)
        │
        ▼  UMAP projection (2D)
        │
        ▼  Wilcoxon rank-sum test → marker genes per cluster
        │
        ▼  Cell type annotation → pbmc3k.h5ad
```

---

## Methods

### Quality Control

![QC Violin](figures/01_qc_violin.png)

Three metrics were assessed per cell:

| Filter | Threshold | Purpose |
|---|---|---|
| Min genes per cell | > 200 | Remove empty droplets |
| Max genes per cell | < 2,500 | Remove likely doublets |
| % mitochondrial reads | < 5% | Remove dying/stressed cells |

> When a cell is dying, cytoplasmic mRNA leaks out while mitochondrial RNA (enclosed in mitochondria) is retained. High MT% reliably flags low-quality or apoptotic cells.

![QC Scatter](figures/02_qc_scatter.png)

---

### Normalization

Library-size normalization was applied by scaling each cell to **10,000 total counts**, followed by log(x+1) transformation to compress dynamic range and stabilize variance across genes.

---

### Highly Variable Gene Selection

![HVG](figures/03_highly_variable.png)

Genes were ranked by **mean expression** and **normalized dispersion** (variance relative to mean). Parameters used:

```python
sc.pp.highly_variable_genes(
    adata,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5
)
```

This reduced the feature space from 32,738 genes → **~1,838 highly informative genes**, dramatically speeding up computation while retaining biological signal.

---

### Dimensionality Reduction

![PCA](figures/04_pca_variance.png)

**PCA** was applied to the scaled HVG matrix with 40 principal components. The variance ratio plot confirms the first 10–15 PCs capture the majority of biological variance.

**UMAP** was then computed from the KNN neighbor graph (k=10, n_pcs=40), projecting cells into 2D space while preserving local neighborhood structure.

---

### Marker Gene Expression on UMAP

![UMAP Markers](figures/05_umap_markers.png)

Three canonical PBMC marker genes overlaid on UMAP:

- **CST3** — monocyte/dendritic cell marker → strongly expressed in the left cluster
- **NKG7** — NK/cytotoxic T cell marker → bottom-right main cluster
- **PPBP** — platelet marker → isolated bottom-right island

---

### Clustering

![UMAP Leiden](figures/06_umap_leiden.png)

The **Leiden algorithm** was applied at `resolution=0.9`, yielding **8 distinct clusters**. Leiden guarantees well-connected communities, improving on the older Louvain approach.

```python
sc.tl.leiden(
    adata,
    resolution=0.9,
    random_state=0,
    flavor="igraph",
    n_iterations=2,
    directed=False,
)
```

---

### Marker Gene Identification & Dotplot

![Dotplot](figures/07_dotplot.png)

Differentially expressed genes were identified per cluster using the **Wilcoxon rank-sum test**. The dotplot shows:

- **Dot size** — fraction of cells in the cluster expressing the gene
- **Color intensity** — mean expression level in the cluster

Each cluster shows distinct, non-overlapping marker patterns — confirming clean separation of cell populations.

---

## Results — Cell Type Annotations

![Annotated UMAP](figures/08_umap_annotated.png)

| Cluster | **Cell Type** | *Key Markers* | Biology |
|---|---|---|---|
| 0 | CD4 T cells | *IL7R, CCR7* | Helper T cells — orchestrate adaptive immune responses |
| 1 | CD14 Monocytes | *LYZ, CD14* | Classical monocytes — phagocytosis & cytokine production |
| 2 | B cells | *CD79A, MS4A1* | Antibody-producing lymphocytes |
| 3 | CD8 T cells | *CD8A, CD8B* | Cytotoxic T cells — kill virally infected cells |
| 4 | NK cells | *GNLY, NKG7* | Natural killer cells — innate cytotoxicity |
| 5 | CD14 Monocytes | *S100A8, LGALS3* | Inflammatory monocyte subpopulation |
| 6 | Dendritic cells | *FCER1A, CST3* | Professional antigen-presenting cells |
| 7 | FCGR3A Monocytes | *FCGR3A, MS4A7* | Non-classical patrolling monocytes |

---

## ⚙️ How to Reproduce

### Install Dependencies

```bash
pip install scanpy anndata matplotlib seaborn pandas numpy scipy
pip install python-igraph leidenalg
```

### Run in Google Colab

1. Upload `notebooks/pbmc3k_analysis.py` to Google Colab
2. Run all cells in order
3. Figures will be saved automatically

### Key Code Snippet

```python
import scanpy as sc

adata = sc.datasets.pbmc3k()
adata.var_names_make_unique()

# QC
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
adata = adata[adata.obs.n_genes_by_counts < 2500]
adata = adata[adata.obs.pct_counts_mt < 5]

# Normalize & cluster
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata, resolution=0.9)
sc.tl.umap(adata)
```

---

## 🔗 References

- [Scanpy Documentation](https://scanpy.readthedocs.io)
- [Scanpy PBMC3k Tutorial](https://scanpy.readthedocs.io/en/latest/tutorials/basics/clustering-2017.html)
- *Wolf et al. (2018)* — Scanpy: large-scale single-cell gene expression data analysis. **Genome Biology**

# GEM-to-Seurat Spatial Transcriptomics Pipeline

This script converts **Stereo-seq GEM files** into **Seurat spatial objects**, merges multiple samples, performs normalization and integration, and runs standard dimensional reduction and clustering.

It is designed for **custom binning–based spatial analysis** when raw Stereo-seq data are available as GEM files but not as native Seurat image objects.

---

## Overview

The script performs the following tasks:

1. **Read Stereo-seq GEM files**
2. **Aggregate reads into spatial bins** (user-defined bin size)
3. **Convert GEM data into Seurat objects**
4. **Construct pseudo-spatial metadata** (no real histology image required)
5. **Merge multiple samples**
6. **Normalize with SCTransform**
7. **Run PCA and UMAP**
8. **Integrate samples using Harmony**
9. **Perform clustering**
10. **Generate QC plots**

This pipeline is suitable for:
- Region-of-interest GEM files
- Downsampled spatial data
- Multi-sample spatial comparisons without registered images

---

## Requirements

### R Packages

The script requires the following R packages:

```r
library(Seurat)
library(data.table)
library(Matrix)
library(jsonlite)
library(ggplot2)

# sample_processing_lasso.py

Interactive ROI selection and GEM/image export for Stereo-seq data using the `stereo` Python API.

## 1. Overview

This script:

1. Reads a tissue gef file (**_tissue.gef**).
2. Runs basic QC on the dataset.
3. Opens an **interactive spatial scatter** view that lets you draw a polygon (lasso) to select a region of interest (ROI).
4. Exports:
   - A **GEM file** containing only the selected ROI.
   - A **cropped TIFF image** corresponding to the selected ROI.

This is useful when you want to downsample a large tissue to one or more manually curated regions for downstream analysis.

## 2. Requirements

- Python 3.8
- The `stereo` package (StereoPy / stomics API), with:
  - `st.io.read_gef`
  - `data.tl.cal_qc`
  - `data.plt.interact_spatial_scatter`
- A working graphical backend (Jupyter notebook, or a Python environment that can open interactive plots).

Typical install (example):

```bash
conda create -n stereo python=3.10
conda activate stereo
pip install stereo

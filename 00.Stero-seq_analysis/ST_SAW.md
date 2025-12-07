# ST_SAW.sh

Stereo-seq SAW v7.0 pipeline script for processing a single spatial transcriptomics sample.

## Overview

`ST_SAW.sh` is a Bash script that runs the **SAW (Stereo-seq Analysis Workflow) v7.0** pipeline for one Stereo-seq sample.  
It performs read mapping, coordinate decoding, gene counting, and image integration to produce standard SAW outputs such as GEF files and QC metrics.

This script is intended as:
- A **reproducible record** of how a sample was processed
- A **template** for running additional samples
- A **methods-level reference** for publications or collaborations

## Requirements

- Linux environment with Bash
- SAW v7.0 (installed locally or via Singularity/Apptainer)
- Adequate CPU and memory resources (HPC recommended)
- Input data and references:
  - Stereo-seq FASTQ files (R1/R2)
  - Mask file from the chip
  - Tissue image (.tar.gz) and image record file (.ipr)
  - Genome reference index (SAW-compatible)
  - Gene annotation GTF file

## Script Structure

The script contains three logical sections:

1. **Sample and path configuration**
2. **Environment / parameter setup**
3. **SAW execution command**

If `set -euo pipefail` is enabled, the script will terminate immediately on errors or missing variables to ensure reproducibility.

## Key Variables

You must edit the following variables before running:

```bash
# ---- Sample-specific settings ------------------------------------------------
SAMPLE_ID=...

# Output directory for this sample
OUT_DIR=...

# Stereo-seq image and mask files
MASK_FILE=...
IMAGE_RECORD_FILE=...
IMAGE_TAR_GZ=...

# ---- Resolve FASTQ file paths -----------------------------------------------
FQ1=...
FQ2=...

# ---- Reference genome and annotation ----------------------------------------
GENOME_INDEX_DIR=...      # STAR index used by SAW
GTF_FILE=...            # gene annotation (GTF)

# ---- SAW / Singularity settings ---------------------------------------------
SAW_SIF=...
THREADS=32
GENOME_SIZE=5          # passed to SAW -genomeSize
SPLIT_COUNT=1             # number of FASTQ splits for parallelization
SPECIES_NAME=...
TISSUE_TYPE=...


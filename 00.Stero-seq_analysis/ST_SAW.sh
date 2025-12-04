#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Stereo-seq analysis with SAW v7.0
# Example command line for one sample, suitable for methods / code sharing
###############################################################################

# ---- Sample-specific settings ------------------------------------------------
SAMPLE_ID="C02647B3"

# Output directory for this sample
OUT_DIR="/path/to/results/${SAMPLE_ID}"

# FASTQ directory (R1/R2 files generated from Stereo-seq sequencing)
FQ_DIR="/path/to/fastq/${SAMPLE_ID}"


# Stereo-seq image and mask files
MASK_FILE="/path/to/mask/${SAMPLE_ID}.barcodeToPos.h5"
IMAGE_RECORD_FILE="/path/to/image/${SAMPLE_ID}.imageRecord.txt"   # example
IMAGE_TAR_GZ="/path/to/image/${SAMPLE_ID}.image.tar.gz"           # example

# ---- Reference genome and annotation ----------------------------------------
GENOME_INDEX_DIR="/path/to/genome_index"      # STAR index used by SAW
GTF_FILE="/path/to/annotation.gtf"            # gene annotation (GTF)

# ---- SAW / Singularity settings ---------------------------------------------
SAW_SIF="/path/to/saw_v7.0.0.sif"
THREADS=32
GENOME_SIZE=5          # passed to SAW -genomeSize
SPLIT_COUNT=1             # number of FASTQ splits for parallelization
SPECIES_NAME="ZH11_T2T"
TISSUE_TYPE="Infl"

# ---- Resolve FASTQ file paths -----------------------------------------------
FQ1="${FQ_DIR}/${SAMPLE_ID}_R1.fastq.gz"
FQ2="${FQ_DIR}/${SAMPLE_ID}_R2.fastq.gz"

# ---- Environment variables for Singularity / HDF5 ---------------------------
export SINGULARITY_BIND="${FQ_DIR},${GENOME_INDEX_DIR},${OUT_DIR},$(dirname "${MASK_FILE}"),$(dirname "${IMAGE_TAR_GZ}")"
export HDF5_USE_FILE_LOCKING=FALSE

# ---- Run SAW pipeline -------------------------------------------------------
sh /path/to/stereoPipeline_v7.0_fromRegister.sh \
  -genomeSize "${GENOME_SIZE_GB}" \
  -splitCount "${SPLIT_COUNT}" \
  -maskFile "${MASK_FILE}" \
  -fq1 "${FQ1}" \
  -fq2 "${FQ2}" \
  -speciesName "${SPECIES_NAME}" \
  -tissueType "${TISSUE_TYPE}" \
  -refIndex "${GENOME_INDEX_DIR}" \
  -annotationFile "${GTF_FILE}" \
  -imageRecordFile "${IMAGE_RECORD_FILE}" \
  -imageCompressedFile "${IMAGE_TAR_GZ}" \
  -sif "${SAW_SIF}" \
  -threads "${THREADS}" \
  -doCellBin N \
  -rRNAremove N \
  -outDir "${OUT_DIR}"
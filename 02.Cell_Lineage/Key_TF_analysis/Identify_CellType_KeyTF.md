# Key TF Spatial Transcriptomics Pipeline

This R program is used to identify the candidate key TF for the newly emerged cell types based on cellular trajectories. The analysis was performed following the protocol described by Qiu et al. with minor adaptations.

---
## Usage
The program performs based on the provided input parameters,as shown below:

Options:
        -i CHARACTER, --input_rds=CHARACTER
                scRNA-seq数据聚类注释后的rds文件,要求细胞类型的列名为'celltype',必须参数.

        -c CHARACTER, --celltype_trajectory=CHARACTER
                细胞分化轨迹关系文件, 主要包含3列: ancestral celltype, new celltype, sister celltype. 必须参数.

        -t CHARACTER, --tf_list=CHARACTER
                转录因子列表文件, 包含1列: TF gene id. 必须参数.

        -l LOGFC_THRESHOLD, --logfc_threshold=LOGFC_THRESHOLD
                FindMarkers()函数中的表达差异倍数logfc.threshold值[default 0.1].

        -w WORKERS, --workers=WORKERS
                Number of workers for parallel processing [default 1];并行计算时使用的核数.

        -o CHARACTER, --output_dir=CHARACTER
                输出结果文件的目录[默认: 当前目录]

        -h, --help
                Show this help message and exit


This pipeline is suitable for:
- Identify the positive and negative candidate key TFs for cell-type specification.

---

## Requirements

### R Packages

The script requires the following R packages:

```r
library(Seurat)
library(dplyr)  ## 管道符
library(optparse)  ## 传入参数解析和帮助文档生成
library(future) ## 多线程任务并行


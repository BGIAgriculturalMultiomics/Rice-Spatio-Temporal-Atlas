#!/usr/bin/env Rscript

## Usage:
##   Rscript go_enrich.R <InPut.list> <gene.fa.iprscan.gene.GO> <Out>

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)   # for barplot()/dotplot()
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript go_enrich.R <InPut.list> <gene.fa.iprscan.gene.GO> <Out>\n")
}

geneListFile <- args[1]   # InPut.list
iprFile      <- args[2]   # gene.fa.iprscan.gene.GO
outPrefix    <- args[3]   # Out

pvalueCutoff <- 0.05

## ------------------------------------------------------------
## 1. Parse iprscan GO file -> <Out>.all.fun.go
## ------------------------------------------------------------

ipr <- read.table(
  iprFile,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

# Expected structure:
# V1 = gene ID
# V2 = number of GO terms (n)
# V3..V(2+n) = "GO:xxxx; description"
res_list <- vector("list", nrow(ipr))
idx <- 1L

for (i in seq_len(nrow(ipr))) {
  id <- ipr[i, 1]
  n  <- ipr[i, 2]

  if (is.na(n) || n <= 0) next

  entries  <- ipr[i, 2 + seq_len(n)]
  go_split <- strsplit(as.character(entries), ";", fixed = TRUE)

  for (g in go_split) {
    if (length(g) < 2) next
    go_id <- g[1]
    desc  <- g[2]

    # clean description (like in the original Perl code)
    desc <- gsub('"', "", desc, fixed = TRUE)
    desc <- gsub("'", "", desc, fixed = TRUE)
    desc <- trimws(desc)

    res_list[[idx]] <- data.frame(
      ID          = id,
      GO          = go_id,
      Description = desc,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
}

res_list <- res_list[!vapply(res_list, is.null, logical(1))]
go_anno_data <- if (length(res_list) > 0) {
  do.call(rbind, res_list)
} else {
  stop("No GO annotations parsed from ipr file.\n")
}

outGoFile <- paste0(outPrefix, ".all.fun.go")
write.table(
  go_anno_data,
  file = outGoFile,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

## ------------------------------------------------------------
## 2. GO enrichment with clusterProfiler::enricher
## ------------------------------------------------------------

# TERM2GENE and TERM2NAME from the annotation we just built
go2gene <- go_anno_data[, c("GO", "ID")]
go2name <- go_anno_data[, c("GO", "Description")]

# gene list
gene_select <- read.table(geneListFile, header = FALSE, stringsAsFactors = FALSE)

# Significant enrichment
ego <- enricher(
  gene_select$V1,
  TERM2GENE    = go2gene,
  TERM2NAME    = go2name,
  pvalueCutoff = pvalueCutoff,
  pAdjustMethod = "fdr",
  qvalueCutoff = 0.2,
  minGSSize    = 0,
  maxGSSize    = 3000
)

# All enrichment
all_res <- enricher(
  gene_select$V1,
  TERM2GENE    = go2gene,
  TERM2NAME    = go2name,
  pvalueCutoff = 1,
  pAdjustMethod = "fdr",
  qvalueCutoff = 1,
  minGSSize    = 0,
  maxGSSize    = 3000
)

# Write result tables
write.table(
  as.data.frame(all_res),
  paste0(outPrefix, ".go_all.xls"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  as.data.frame(ego),
  paste0(outPrefix, ".go_enrichment.significant.xls"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

## ------------------------------------------------------------
## 3. Plots
## ------------------------------------------------------------

pdf(paste0(outPrefix, ".go_enrichment.all.V1.pdf"), width = 10, height = 8)
barplot(all_res, showCategory = 20, label_format = 60, color = "p.adjust")
dev.off()

pdf(paste0(outPrefix, ".go_enrichment.all.V2.pdf"), width = 10, height = 8)
dotplot(all_res, showCategory = 20, label_format = 60, color = "p.adjust")
dev.off()

pdf(paste0(outPrefix, ".go_enrichment.sign.V1.pdf"), width = 10, height = 8)
barplot(ego, showCategory = 20, label_format = 60, color = "p.adjust")
dev.off()

pdf(paste0(outPrefix, ".go_enrichment.sign.V2.pdf"), width = 10, height = 8)
dotplot(ego, showCategory = 20)
dev.off()
# GO Enrichment from IPRscan Annotation

`go_enrich.R` performs **Gene Ontology (GO) enrichment analysis** using a custom IPScan-based GO annotation file and a list of genes of interest. It uses **clusterProfiler** for enrichment and **enrichplot** for visualization.

---

## 1. Usage

```bash
Rscript go_enrich.R <InPut.list> <gene.fa.iprscan.gene.GO> <Out>

library(Seurat)
library(sceasy)
library(data.table)
library(readxl)
library(tidyverse)
library(anndata)
library(SeuratDisk)
library(PlantPhoneDB)
library(UniprotR)
library(riceidconverter)
library(dplyr)

library(ggplotify)
library(circlize)
library(RColorBrewer)
library(igraph)
library(pheatmap)
library(plyr)

data = readRDS("./data.rds")

Idents(data) = 'celltype_L2'

data <- SCTransform(data,assay = "RNA",verbose = FALSE)

# saveRDS(data, './sct.rds')

LR_pair_anno = read.csv('./from_ipSAE/H025_5000/LRs_ipSAE.xls', sep='\t')
LR_pair = LR_pair_anno[, c('Ligands', 'Receptors')]
print(dim(LR_pair))
# head(LR_pair, 3)

LRp <- LRscore(data@assays$SCT@data, LRdb=LR_pair, cluster = Idents(data), min.pct = 0.01,iterations=100, method='Average')

outpath = './results/'
write.csv(LRp, file=paste0(outpath, "LRscore.csv"), quote=F, row.names=F)

LRp_sig <- LRp %>%
    filter(Pvalue<0.05) 


write.csv(LRp_sig, file=paste0(outpath, "LRscore_sig.csv"), quote=F, row.names=F)

df.net = read.csv(paste0(outpath, 'LRscore_sig.csv'))
head(df.net, 3)

source_counts <- table(df.net$Ligands_cell)
target_counts <- table(df.net$Receptors_cell)
# 将两个向量合并为一个数据框
data <- data.frame(
  CellType = unique(c(names(source_counts), names(target_counts)))
  # Count = c(source_counts, target_counts)
  # Category = factor(rep(c("Source", "Target"), each = length(target_counts)))
)
data$Source <- ifelse(data$CellType %in% names(source_counts), source_counts[data$CellType], 0)
data$Target <- ifelse(data$CellType %in% names(target_counts), target_counts[data$CellType], 0)

long_df <- data %>%
  pivot_longer(
    cols = c(Source, Target),
    names_to = "Category",
    values_to = "Count"
  )

p = ggplot(long_df, aes(x = CellType, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Cell Type Counts by Source and Target", x = "Cell Type", y = "Count") + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(p, width = 10, height = 6, filename = paste0(outpath, 'lrs_stat_bar.pdf'))

out_dir = outpath

interaction_count <- LRp_sig %>%
    group_by(Ligands_cell, Receptors_cell) %>%
    dplyr::summarise(Number=n(), .groups = 'drop')

interaction_count %>%
    mutate(Type=ifelse(Ligands_cell==Receptors_cell,"Autocrine","Paracrine")) %>%
    group_by(Type) %>%
    summarise(Number=sum(Number))

Autocrine <- interaction_count[interaction_count$Ligands_cell==interaction_count$Receptors_cell,]
Paracrine <- interaction_count[interaction_count$Ligands_cell !=interaction_count$Receptors_cell,]

cs = unique(append(interaction_count$Ligands_cell, interaction_count$Receptors_cell))
length(cs)

mycolor <- c('#d60000','#018700','#b500ff','#05acc6','#97ff00','#ffa52f','#ff8ec8','#79525e',
                '#00fdcf','#afa5ff','#93ac83','#9a6900','#366962','#d3008c','#fdf490','#c86e66',
              '#9ee2ff','#00c846','#a877ac','#b8ba01','#f4bfb1', "tan1","darkred","salmon","#6B8E23",'#F7F7F7', 
            '#d60000','#018700','#b500ff','#05acc6','#97ff00','#ffa52f','#ff8ec8','#79525e','#00fdcf')

# options(repr.plot.width=14, repr.plot.height=14)
pdf(paste0(out_dir, "/CCI_circle.pdf", sep=''), width=14, height=14)
CCI_circle(Paracrine, mycolor[1:length(cs)])
dev.off()

myCCI_network <- function(interaction_count,mycolor, vertex.label.cex=1.5, edge.label.cex=1, title="",edgeLabel=TRUE){
    color <- data.frame(
        Ligands_cell=unique(interaction_count$Ligands_cell),
        color=mycolor[1:length(unique(interaction_count$Ligands_cell))]
    )
    interaction_count <- interaction_count %>%
        inner_join(color)
    net <- graph_from_data_frame(interaction_count)
    karate_groups <- cluster_optimal(net)
    coords <- layout_in_circle(net, order= order(membership(karate_groups))) 
    E(net)$width  <- E(net)$Number/10 
    V(net)$color <- mycolor[get.vertex.attribute(net)$name]
    #E(net)$color <- mycolor
    if(edgeLabel){
        E(net)$label <- E(net)$Number
    }
    pic <- plot(net, edge.arrow.size=1, 
     edge.curved=0.2,
     vertex.label.color="black",
     layout = coords,
  edge.label.cex= edge.label.cex,
     vertex.label.cex=vertex.label.cex,main=title)
}

options(repr.plot.width=10, repr.plot.height=10)
pdf(paste0(out_dir, "/CCI_network.pdf", sep=''), width=14, height=14)
myCCI_network(interaction_count, mycolor, vertex.label.cex=2, edge.label.cex=1)
dev.off()

options(repr.plot.width=30, repr.plot.height=32)
pdf(paste0(out_dir, "/heatmap_count.pdf", sep=''))
heatmap_count(interaction_count,text_size=10,number_size=2,decimal=4,title="Heatmap")
dev.off()

Top10 <- LRp_sig %>%
    arrange(desc(Score)) %>%
    select(LR_pair) %>%
    unique() %>%
    head(10) %>%
    inner_join(LRp) %>%
    select(LR_pair,Cell_pair,Score) %>%
    unique() %>%  ### 自己加的
    spread(.,Cell_pair,Score) %>%
    replace(is.na(.), 0)
# Top10
rownames(Top10) <- Top10$LR_pair
Top10 <- Top10[,-1]
Top10 <- t(Top10)

# options(repr.plot.width=8, repr.plot.height=25)
# pdf("./05.PlantPhoneDB/output/Top10_heatmap.pdf")
print(length(Top10))
pic6 <- pheatmap(
    Top10, scale="none",angle_col=90,fontsize_row=8,
    cluster_rows = T,cluster_cols = F,show_colnames=T, 
    width=6, height=12,
    filename = paste0(out_dir, "/Top10_heatmap_symbol.pdf", sep='')
    # filename = paste0(out_dir, "/Top10_heatmap.pdf", sep='')
)

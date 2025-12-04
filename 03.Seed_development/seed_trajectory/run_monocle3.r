suppressPackageStartupMessages({
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(dplyr)
library(data.table)
library(Matrix)
library(rjson)
library(RColorBrewer)
library(ggplot2)
    library(Seurat)
    library(dplyr)
    library(SingleR)
    library(pheatmap)
    library(ggplot2)
    library(RColorBrewer)
    library(future)
    library(dplyr)
    library(patchwork)
    library(tidyverse)
    library("harmony")
    library(monocle3)
})

set.seed(123)


source("/share/appspace_data/shared_groups/BGI/04.Project/hewmRice/90.Bin/03.Fun/aa.r")


readRDS("05.Day05_Day15.rds")->seob


cluster_Palette2 <- c('dodgerblue2', '#E31A1C', 'green4', '#6A3D9A', '#FF7F00', 'black', 'gold1',
                         'skyblue2', '#FB9A99', 'palegreen2', '#CAB2D6', '#FDBF6F', 'gray70', 'khaki2',
                         'maroon', 'orchid1', 'deeppink1', 'blue1', 'steelblue4', 'darkturquoise',
                         'green1','yellow4', 'yellow3','darkorange4',
                      'steelblue2','darkolivegreen3','darkorchid3','lightcoral','mediumpurple2')

iPlotAA <- function(object, features, pt.size = 0.02){
    plot <- ggplot(object, aes_string(x = 'newx', y = 'newy', color = features)) +
            geom_point(shape = 19, size = pt.size) + theme_void() +coord_fixed()+
            theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
#axis.text: ?𹰿𽳞杞存𻓞??? axis.ticks: ?𹰿𽳞杞村𷰿搴?嚎锛𼠠 panel.grid: 缃𺀿𽰿绾?
                  axis.title = element_blank(), axis.line = element_blank(), plot.margin=margin(t=1,b=1,r=1,l=1))
     plot <- plot + scale_color_manual(values = cluster_Palette2) +
                guides(colour = guide_legend(override.aes = list(size=3), nrow = 10))
    
    plot <- plot + theme_void()
    return(plot)
}

iPlotCC <- function(object, features, pt.size = 0.88){
    plot <- ggplot(object, aes_string(x = 'newxCC', y = 'newyCC', color = features)) +
            geom_point(shape = 19, size = pt.size) + theme_void() +coord_fixed()+
            theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
#axis.text: ?𹰿𽳞杞存𻓞??? axis.ticks: ?𹰿𽳞杞村𷰿搴?嚎锛𼠠 panel.grid: 缃𺀿𽰿绾?
                  axis.title = element_blank(), axis.line = element_blank(), plot.margin=margin(t=1,b=1,r=1,l=1))
     plot <- plot + scale_color_manual(values = cluster_Palette2) +
                guides(colour = guide_legend(override.aes = list(size=3), nrow = 10))
    
    plot <- plot + theme_void()
    return(plot)
}


head(seob,2)


unique(seob@meta.data$Day)
seob@meta.data->df
#write.table(file="01ALL.table",df,sep="\t", quote=FALSE, row.names=TRUE)
options(repr.plot.width =14,repr.plot.height =9)
iPlotCC(df,'ClusterDD')


DayyNow <- subset(seob, subset = Day %in% c("Day08"))


iPlotCC(DayyNow@meta.data,'ALLAnoV1',pt.size = 0.82)


iPlotCC(DayyNow@meta.data,'ClusterDD',pt.size = 0.82)


head(DayyNow,2)


DayyNowPeiRu <- subset(DayyNow, subset = ClusterDD  %in% c("DD01",'DD02',"DD03","DD04","DD05","DD07"))

   #6736 Peripheral-endosperm_05-08DAP_1	DD01
   #4326 Subaleurone_05-15DAP	DD03
   #3454 Endosperm_05-15DAP	DD04
   #2332 Peripheral-endosperm_05-08DAP_2	DD07
   #4874 Starchy-endosperm_05-15DAP	DD02


   iPlotCC(DayyNowPeiRu@meta.data,'ALLAnoV1',pt.size = 0.88)


   heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
library(viridis)
#scale_color_viridis(option = "C", begin = 0,end = 0.9, na.value = '#f2eee7') 

iPlotV2 <- function(seob, color_attr, pt.size = 0.2){
    objectA <- data.frame('x' = seob[[]]$x, 'y' = seob[[]]$y,  'ColorAttr' = seob[[]]$ColorAttr, 'color_attr' = seob[[]][[color_attr]])
    #print(head(objectA,3))
    object  <- objectA[objectA$color_attr >= 0, ]
    plot <- ggplot(object, aes_string(x = 'x', y = 'y', color = 'color_attr')) +
            geom_point(shape = 19, size = pt.size) + theme_void() +coord_fixed()+
            theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
                  axis.title = element_blank(), axis.line = element_blank(), plot.margin=margin(t=1,b=1,r=1,l=1))
    if (color_attr %in% c('pseudotime','PeiRuPseudotime','PeiPseudotime')){
        plot <- plot + scale_color_gradientn(colours = heatmap_Palette(100)) +
        guides(colour = "colorbar", color=guide_legend(override.aes = list(size=1), nrow = 10, title = "pseudotime")) 
    }else if(color_attr %in% c('seurat_clusters')){
        plot <- plot + scale_color_manual(values = cluster_Palette) +
                guides(colour = guide_legend(override.aes = list(size=3), nrow = 10))
    }
    plot <- plot + theme_void()
    return(plot)
}

#iPlotV2(seob, color_attr = 'PeiRuPseudotime', pt.size = 0.05)


visualALL <- function(seob) {
  df <- data.frame('x' = seob[[]]$x, 'y' = seob[[]]$y, 'cluster' = seob[[]]$seurat_clusters,'sample'= seob[[]]$sample )
 # colorPalette <- ColorPalette(16)
  point_size = 0.1
  p0 <-
    ggplot(df, aes(x = x, y = y, color = cluster)) +
    geom_point(shape = 19, size = point_size) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank(),
      legend.position = 'right'
    ) +
    scale_color_manual(values = cluster_Palette) +
    guides(colour = guide_legend(
      override.aes = list(size = 3),
      nrow = 15,
      title = 'Clusters'
    ))  +
    theme_void()+coord_fixed()

  p0 <- p0 + facet_grid(~sample)
# p0 <- p0 + theme(aspect.ratio = 1,      plot.margin = margin(10, 10, 10, 10)) +coord_fixed()
  
  return (p0)
}
#visualALL(PeiBB)


VisualGeneExpSTOV2 <- function(MakerGeneList,  seob)
{
    gene_expression <- FetchData(object = seob, vars = MakerGeneList, slot = "data")
    Tmpseob <- AddMetaData(object = seob, metadata = gene_expression)
    plots <- list()
    for (i in 1:length(MakerGeneList)){
        data <- data.frame(Tmpseob@meta.data)
        plot_col <- c('x','y',MakerGeneList[i],'sample')
        df <- data[,colnames(data) %in% plot_col]
        colnames(df) <- c('x','y','Day','exp')
        point_size = 0.3
        plots[[i]] <- ggplot(df, aes(x = x, y = y, color = exp)) +
        geom_point(shape = 19, size = point_size) +  ## you can modify the point size to avoid overlap
        theme(
              #      axis.text = element_blank(),
              #      axis.ticks = element_blank(),
              panel.grid = element_blank(),
              axis.title = element_blank(),
              axis.line = element_blank(),
              legend.position = 'right'
              ) +
coord_fixed()  +
#scale_fill_gradientn(colours = colorPalette)+
# '#0000ff','#ffffff','#ff0000'
scale_colour_gradient2(low = "#0000ff", mid = "#ffffff", high = "#ff0000") + # ??
guides(colour = guide_legend(
                             override.aes = list(size = 3),
                             nrow = 15,
                             title = 'expression'
                             ))  +
theme_void()+ggtitle(MakerGeneList[i]) +facet_grid(~ Day)
    }
    return (plots)
}

#TargetMakerGene<-c("OsZH11G0100065600.01")  # "OsZH11G0818334700.01"
#VisualGeneExpSTOV2(TargetMakerGene,dayall)


saveRDS(DayyNowPeiRu, file = "10.Day08PeiRuV2.rds")


DayyNowBB<- subset(seob, subset = Day %in% c("Day15"))


iPlotCC(DayyNowBB@meta.data,'ALLAnoV1',pt.size = 0.82)
iPlotCC(DayyNowBB@meta.data,'ClusterDD',pt.size = 0.82)


DayyNowPeiRuBB  <- subset(DayyNowBB, subset = ClusterDD  %in% c("DD01",'DD02',"DD03","DD05","DD06"))
iPlotCC(DayyNowPeiRuBB@meta.data,'ALLAnoV1',pt.size = 0.88)


saveRDS(DayyNowPeiRuBB, file = "10.Day15PeiRuV2.rds")


readRDS("10.Day08PeiRuV2.rds")->DayyNowPeiRu


DayyNowPeiRuBB  <- subset(DayyNowPeiRu, subset = ClusterDD  %in% c("DD04","DD05","DD02"))
iPlotCC(DayyNowPeiRuBB@meta.data,'ALLAnoV1',pt.size = 0.88)


sce<-DayyNowPeiRuBB
#创建CDS对象并预处理数据
data <- GetAssayData(sce, assay = 'RNA', slot = 'counts')
cell_metadata <- sce@meta.data
head(cell_metadata,3)
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
# 创建CDS(Cell Data Set)对象
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
   cds <- preprocess_cds(cds, num_dim = 50)
# 要对批次进行校正,我们这是整合  须要处理批次
#cds = preprocess_cds(cds, num_dim = 100, residual_model_formula_str="~ sample")
cds <- align_cds(cds,  preprocess_method = c("PCA"),  alignment_group ="MSample",  alignment_k = 50,  residual_model_formula_str = "~ MSample")
#cds <- align_cds(cds,  preprocess_method = c("PCA"),  alignment_group ="~sample",  alignment_k = 50)
head(cds,3)


options(repr.plot.width =8,repr.plot.height =8)
plot_pc_variance_explained(cds)


#umap降维
cds <- reduce_dimension(cds, preprocess_method = "Aligned")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="ALLAnoV1", group_label_size = 8) + ggtitle('cds.umap')
##从seurat导入整合过的umap坐标
p1


cds.embed <- cds@int_colData$reducedDims$UMAP


#int.embed <- Embeddings(sce, reduction = "UMAP_1")
#int.embed <- int.embed[rownames(cds.embed),]
#cds@int_colData$reducedDims$UMAP <- int.embed
#p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="ALLAnoV1", group_label_size = 8) + ggtitle('int.umap')
#options(repr.plot.width =11,repr.plot.height =6)
#p1+p2


cds <- cluster_cells(cds)
#saveRDS(cds,"02.PeiRuCDS.rds")


cds <- learn_graph(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE,group_label_size=8,color_cells_by="ALLAnoV1") + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "MSample", show_trajectory_graph = FALSE,group_label_size=8) + 
        ggtitle("label by sample")
p = wrap_plots(p1, p2)
p


p = plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, ,group_label_size=8,color_cells_by="ALLAnoV1",
               label_branch_points = FALSE)+scale_color_manual(values = cluster_Palette) +
                guides(colour = guide_legend(override.aes = list(size=3), nrow = 10))
#+scale_colour_gradient(colours = rev(brewer.pal(n = 11, name = "Spectral")))
p


get_earliest_principal_node <- function(cds, time_bin="EAS_05-15DAP")
{ 
  cell_ids <- which(colData(cds)[, "ALLAnoV1"] == time_bin)
    
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

#get_earliest_principal_node(cds)

cds <- order_cells(cds, root_pr_nodes=c(get_earliest_principal_node(cds,"EAS_05-15DAP")))
plot_cells(cds,color_cells_by = "pseudotime",label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5,group_label_size=8)+scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))


pdf("Day08UmapTime.pdf")
plot_cells(cds,color_cells_by = "pseudotime",label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5,group_label_size=8)+scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) 
dev.off()
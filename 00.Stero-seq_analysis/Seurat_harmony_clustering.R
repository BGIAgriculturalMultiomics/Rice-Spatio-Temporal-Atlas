library(Seurat)

gem_file_sample1 = /path/to/gem_file_sample1
gem_file_sample2 = /path/to/gem_file_sample2
gem_file_sample3 = /path/to/gem_file_sample3

gem_to_seuratObject <- function(data, prefix = 'sample', binsize = 100){
  #' read the gem file  
  data <- fread(file = data)
  data$x <- as.numeric(data$x)
  data$y <- as.numeric(data$y)
  
  #' group counts into bins
  data$x <- trunc(data$x / binsize) * binsize
  data$y <- trunc(data$y / binsize) * binsize
  
  if ('MIDCounts' %in% colnames(data)) {
    data$MIDCounts <- as.numeric(data$MIDCounts)
    data <- data[, .(counts=sum(MIDCounts)), by = .(geneID, x, y)]
  } else if ('UMICount' %in% colnames(data)) {
    data$UMICount <- as.numeric(data$UMICount)
    data <- data[, .(counts=sum(UMICount)), by = .(geneID, x, y)]
  } else if ('MIDCount' %in% colnames(data)) {
    data$MIDCount <- as.numeric(data$MIDCount)
    data <- data[, .(counts=sum(MIDCount)), by = .(geneID, x, y)]
  }
  
  #' create sparse matrix from stereo
  data$cell <- paste0(prefix, ':', data$x, '-', data$y)
  data$geneIdx <- match(data$geneID, unique(data$geneID))
  data$cellIdx <- match(data$cell, unique(data$cell))
  
  mat <- sparseMatrix(i = data$geneIdx, j = data$cellIdx, x = data$counts,
                      dimnames = list(unique(data$geneID), unique(data$cell)))
  
  cell_coords <- unique(data[, c('cell', 'x', 'y')])
  
  rownames(cell_coords) <- cell_coords$cell
  
  seurat_spatialObj <- CreateSeuratObject(counts = mat, project = 'Stereo', assay = 'Spatial',
                                          names.delim = ':', meta.data = cell_coords)
  
  
  #' create pseudo image
  cell_coords$x <- cell_coords$x - min(cell_coords$x) + 1
  cell_coords$y <- cell_coords$y - min(cell_coords$y) + 1
  
  tissue_lowres_image <- matrix(1, max(cell_coords$y), max(cell_coords$x))
  
  # tissue_positions_list <- data.frame(row.names = cell_coords$cell,
  #                                     tissue = 1,
  #                                     row = cell_coords$y, col = cell_coords$x,
  #                                     imagerow = cell_coords$y, imagecol = cell_coords$x)
  
  tissue_positions_list <- data.frame(row.names = cell_coords$cell,
                                      tissue = 1,
                                      row = cell_coords$x, col = cell_coords$y,
                                      imagerow = cell_coords$x, imagecol = cell_coords$y)
  
  scalefactors_json <- toJSON(list(fiducial_diameter_fullres = binsize,
                                   tissue_hires_scalef = 1,
                                   tissue_lowres_scalef = 1))

sce.list=list()
sample_names=c("sample_1","sample_2","sample_3")
file_names=c(gem_file_sample1, gem_file_sample2, gem_file_sample3)
for (n in 1:length(sample_names)) {
	obj = gem_to_seuratObject(file_names[n], binsize=50)
	obj$sample = sample_names[n]
	obj$project = "Inflourescence Formation"
	obj$stage = "PBM"
	sce.list = append(sce.list, list(obj))
}

sce=merge(sce.list[[1]],sce.list[2:3])
sce <- SCTransform(sce)
sce <- RunPCA(sce)
sce <- RunUMAP(sce,reduction="pca", reduction.name="umap.unintegrated", n.neighbors=k_snn, dims=1:ndim)
sce <- IntegrateLayers(object = sce, method = HarmonyIntegration, normalization.method = "SCT", verbose = F,orig.reduction="pca",new.reduction="harmony")

#QC
p2=ElbowPlot(sce, ndims = 50, reduction = "pca")
figure_name2=paste0("QC_",sub_sample,"_Elbow.png")
ggsave(figure_name2,p2,bg="white")

#Dimension reduction and clustering
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:30, k.param=15)
sce <- FindClusters(sce,resolution=resolution,graph.name="SCT_snn",cluster.name="harmony_clusters")
sce <- RunTSNE(sce,reduction="harmony",dims=1:30)
sce <- RunUMAP(sce,reduction="harmony", reduction.name="umap", n.neighbors=15, dims=1:30)





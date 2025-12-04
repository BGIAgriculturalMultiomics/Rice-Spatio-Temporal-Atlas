# 插件运行所需的R包
suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)  ## 管道符
    library(optparse)  ## 传入参数解析和帮助文档生成
    library(future) ## 多线程任务并行
})


# 定义程序的简单介绍
program_description <- "
This R program is used to identify key TF condidates for the newly emerged cell types based on cellular trajectories.
Display it along with command-line options when --help is used.

The program performs based on the provided input parameters,as shown below.
"

# 定义命令行参数
option_list <- list(
  make_option(c("-i", "--input_rds"), type="character", default=NULL, 
              help="scRNA-seq数据聚类注释后的rds文件,要求细胞类型的列名为'celltype'.", metavar="character"),
  make_option(c("-c", "--celltype_trajectory"), type="character", default=NULL, 
              help="细胞分化轨迹关系文件, 主要包含3列: ancestral celltype, new celltype, sister celltype.", metavar="character"),
  make_option(c("-t", "--tf_list"), type="character", default=NULL,
	     help="转录因子列表文件, 包含1列: TF gene id.", metavar="character"),
  make_option(c("-l", "--logfc_threshold"), type="numeric", default=0.1,
             help="FindMarkers()函数中的表达差异倍数logfc.threshold值[default %default]."),
  make_option(c("-w", "--workers"), type="integer", default=1, 
              help="Number of workers for parallel processing [default %default];并行计算时使用的核数."),
  make_option(c("-o", "--output_dir"), type="character", default=getwd(), 
              help="输出结果文件的目录[默认: 当前目录]", metavar="character")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list=option_list, description=program_description)
opt <- parse_args(opt_parser)

output_dir <- opt$output_dir

# 检查是否提供了必要的参数
if (is.null(opt$input_rds) || is.null(opt$celltype_trajectory) || is.null(opt$tf_list)) {
  print_help(opt_parser)
  stop("缺少输入文件，必须提供所有参数", call.=FALSE)
}

## 建输出目录 
dir.create("markers")
dir.create("output")

## 读入rds文件
sce <- readRDS(opt$input_rds)

temp <- table(sce$celltype)
print(temp)

### 读取细胞类型间的关系,即new cell type与pseudo-ancestral cell, sisiter cell type.
link <- read.table(opt$celltype_trajectory, sep="\t",header=T)
colnames(link) <- c("Ancestral_cell","New_celltype","Sister_celltype")
link <- link %>% mutate(flag = paste0("N",1:n(),"_vs_A",1:n()))

### ——step 1: 预测celltype间的markers;
### new cell type与假祖先细胞的差异基因：
Idents(sce) <- "celltype"

options(future.globals.maxSize=30000000000)
plan("multisession", workers = opt$workers)
for(i in 1:nrow(link)){
    all.markers <- FindMarkers(sce,
                          assay = 'RNA',
                          slot = 'data',
                          ident.1 = link$New_celltype[i],
                          ident.2 = link$Ancestral_cell[i],
                          logfc.threshold = opt$logfc_threshold,
                          min.pct = 0.1,  # 在其中一组中细胞占比>10%即可.         
                          #only.pos = T
			  )
   #avg_log2FC表示的是基因在ident.1细胞群体中的表达水平相对于ident.2细胞群体的对数2倍变化.
   all.markers.filter <- all.markers[all.markers$p_val_adj<0.05, ]  ## —— 正调控条件1,负调控条件2
   ## 筛选New_celltype中正调控基因(logfc >0) & 本newly emerged cell type的至少10%的细胞中表达(pct.1 > 0.1)
   positive <- all.markers.filter[all.markers.filter$avg_log2FC >0 & all.markers.filter$pct.1 > 0.1, ]  ## ——正调控条件1,3

   ## 筛选New_celltype中负调控基因(logfc <0). pct.2(祖先细胞), pct.2和pct.1不一定都>0.1, 所以加过滤条件pct.1 >0.1 
   negative <- all.markers.filter[all.markers.filter$avg_log2FC <0 & all.markers.filter$pct.2 > 0.1, ] ## ——负调控条件1,2
   write.table(positive, file = paste0(output_dir,"/markers/N",i,"_vs_A",i,".cluster.AllmarkerGene.positive.txt"), sep = "\t")
   write.table(negative, file = paste0(output_dir,"/markers/N",i,"_vs_A",i,".cluster.AllmarkerGene.negative.txt"), sep = "\t")
}

### new cell type与sister cell type的差异基因:
options(future.globals.maxSize=30000000000)
plan("multisession", workers = opt$workers)

for(i in 1:nrow(link)){
    if(is.na(link$Sister_celltype[i]) || link$Sister_celltype[i] == ""){
        message("Sister_celltype of row", i, " is null.") ## 有的new cell type没有sister cell type.
    }else{
    sister_cell_types <- strsplit(link$Sister_celltype[i], ",")[[1]]
    all.markers <- FindMarkers(sce,
                          assay = 'RNA',
                          slot = 'data',
                          ident.1 = link$New_celltype[i],
                          ident.2 = sister_cell_types,
                          logfc.threshold = opt$logfc_threshold,
                          min.pct = 0.01) #在其中一组中细胞占比>1%,先全部输出,表达占比10%后面通过pct.1和pct.2过滤筛选.
                          # only.pos = T) 

   all.markers.filter <- all.markers[all.markers$p_val_adj<0.05, ] ## ——正调控条件2
   
   ## 筛选在New_celltype中正调控的基因 —— for正调控TF条件(2).
   positive <- all.markers.filter[all.markers.filter$avg_log2FC >0, ]  ## ——正调控条件2
   ## 筛选New_celltype中负调控基因(logfc <0). pct.2(siste cell), pct.1不一定>0.1, 所以加过滤条件pct.1 >0.1
   negative <- all.markers.filter[all.markers.filter$avg_log2FC <0 & all.markers.filter$pct.1 > 0.1 & all.markers.filter$pct.2 >0.1,] ##——负调控条件3,4
   write.table(positive, file = paste0(output_dir,"/markers/N",i,"_vs_S",i,".cluster.AllmarkerGene.positive.txt"), sep = "\t")
   write.table(negative, file = paste0(output_dir,"/markers/N",i,"_vs_S",i,".cluster.AllmarkerGene.negative.txt"), sep = "\t")
   }
}


### ——step2 : 预测celltype的正调控key TF和计算score
### 读取TF list
tf <- read.table(opt$tf_list)
colnames(tf) = c("TF")

TF <- unique(tf$TF)
print(length(TF))

### 正调控key TF
#output <- c("Order","Markers_num","keyTF_num","Key_TFs")  ## output title
output1 <- data.frame()
for(i in 1:nrow(link)){
    if(is.na(link$Sister_celltype[i]) || link$Sister_celltype[i] == ""){
        #message("sister_cell of the ", i, " row is null.") ## 有的new cell type没有sister cell type.
    	A <- read.table(paste0(output_dir,"/markers/N",i,"_vs_A",i,".cluster.AllmarkerGene.positive.txt"),sep="\t",header=T)
    	#S <- read.table(paste0("./markers/N",i,"_vs_S",i,".cluster.AllmarkerGene.txt"),sep="\t",header=T)
    	A_genes <- rownames(A)
    	#S_genes <- rownames(S)
    	#key_markers <- intersect(A_genes, S_genes)
    	len1 <- length(A_genes)
    	key_TF <- intersect(A_genes,TF)
    	keyTF_len <- length(key_TF)
    	print(keyTF_len)

    	prefix = paste0("N",i,"_vs_A",i)
    	key_TF <- paste(key_TF,collapse = ",")
    	result <- data.frame(Flag=prefix,Sister_celltype="No",Markers_num=len1,keyTF_num=keyTF_len,TF_IDs=key_TF)
        output1 <- rbind(output1,result)
    }else{
        ancester <- read.table(paste0(output_dir,"/markers/N",i,"_vs_A",i,".cluster.AllmarkerGene.positive.txt"),sep="\t",header=T)
        sister <- read.table(paste0(output_dir,"/markers/N",i,"_vs_S",i,".cluster.AllmarkerGene.positive.txt"),sep="\t",header=T)
    #S <- read.table(paste0("./markers/N",i,"_vs_S",i,".cluster.AllmarkerGene.txt"),sep="\t",header=T)
        ancester_genes <- rownames(ancester)
        sister_genes <- rownames(sister)
        key_markers <- intersect(ancester_genes, sister_genes)
        print(length(ancester_genes))
        print(length(sister_genes))
        len1 <- length(key_markers)
        print(len1)

        key_TF <- intersect(key_markers,TF)
        keyTF_len <- length(key_TF)
        print(keyTF_len)

        prefix = paste0("N",i,"_vs_A",i)
        key_TF <- paste(key_TF,collapse = ",")
        result <- data.frame(Flag=prefix,Sister_celltype="Yes",Markers_num=len1,keyTF_num=keyTF_len,TF_IDs=key_TF)    
        output1 <- rbind(output1,result)
    }
}
write.table(output1,file = paste0(output_dir,"/output/1.positive_NewCelltype_key_TFs_output.txt"), row.names = FALSE, sep="\t")

### 计算正调控key TF的得分
results <- data.frame()
for (i in 1:nrow(output1)){
    if(output1$Sister_celltype[i] == "No"){  ## 有些new cell type是没有sister cell type的
        markers <- read.table(paste0(output_dir,"/markers/N",i,"_vs_A",i,".cluster.AllmarkerGene.positive.txt"),sep="\t",header=T)
        flag_tmp <- paste0("N",i,"_vs_A",i)
        markers$GeneID <- rownames(markers)
        TF_IDs <- output1[output1$Flag == flag_tmp, "TF_IDs"]
        TF_genes <- unlist(strsplit(TF_IDs, ","))
        #print(length(TF_genes))
        if(length(TF_genes)>0){  ## 有的new cell type没有key TFs.
            A_output <- markers[markers$GeneID %in% TF_genes, c("avg_log2FC","pct.1","pct.2","p_val_adj","GeneID")]
            temp_df <- as.data.frame(matrix(NA, nrow = nrow(A_output), ncol = ncol(A_output)))  # 生成一个空的data.frame，用于跟有sister的情况合并；
            colnames(temp_df) <- paste("S", colnames(A_output), sep = ".")
            output <- cbind(A_output, temp_df)
            output$score <- scale(output$avg_log2FC)
            output$flag <- flag_tmp
            output$Ancestral_cell <- link[link$flag == flag_tmp,"Ancestral_cell"] # 添加假祖先cell type
            output$New_celltype <- link[link$flag == flag_tmp,"New_celltype"]  #添加new cell type
            output$Sister_celltype <- "No"
            sorted_df <- output[order(-output$score), ]
            #new_df <- select(sorted_df,"flag","ancestral_cell","New_celltype","GeneID","score","avg_log2FC","pct.1","pct.2","p_val_adj")
            new_df <- select(sorted_df,"flag","Ancestral_cell","New_celltype","Sister_celltype","GeneID","score","avg_log2FC","pct.1","pct.2","p_val_adj","S.avg_log2FC","S.pct.1","S.pct.2","S.p_val_adj")
            names(new_df)[names(new_df)=="GeneID"] <- "key_TF"
            results <- rbind(results,new_df)
        }
    }else{
        print(i)
        flag_tmp <- paste0("N",i,"_vs_A",i)
        ancester_markers <- read.table(paste0(output_dir,"/markers/N",i,"_vs_A",i,".cluster.AllmarkerGene.positive.txt"),sep="\t",header=T)
        sister_markers <- read.table(paste0(output_dir,"/markers/N",i,"_vs_S",i,".cluster.AllmarkerGene.positive.txt"),sep="\t",header=T)
        ancester_markers$GeneID <- rownames(ancester_markers)
        sister_markers$GeneID <- rownames(sister_markers)
        TF_IDs <- output1[output1$Flag == flag_tmp, "TF_IDs"]
        TF_genes <- unlist(strsplit(TF_IDs, ","))
        print(length(TF_genes))
        if(length(TF_genes)>0){  ## 有的new cell type没有key TFs.
            A_output <- ancester_markers[ancester_markers$GeneID %in% TF_genes, c("avg_log2FC","pct.1","pct.2","p_val_adj","GeneID")]
            S_output <- sister_markers[sister_markers$GeneID %in% TF_genes, c("avg_log2FC","pct.1","pct.2","p_val_adj","GeneID")]
            scaled_log_fold_change1 <- scale(A_output$avg_log2FC)  ## 对祖先avg_logFC进行scale
            scaled_log_fold_change2 <- scale(S_output$avg_log2FC)  ## 对sister细胞的avg_logFC进行scale
            matrix_data <- cbind(scaled_log_fold_change1, scaled_log_fold_change2)  # 合并2列;
            row_avg <- apply(matrix_data, 1, function(x) mean(x))  ## 求2个值的平均值。
            colnames(S_output) <- paste("S", colnames(S_output), sep = ".")  # 由于列名相同，在sister markers的列名上加上前缀“S.”
            output <- cbind(A_output, S_output)  ## 把祖先和sister的差异倍数等信息合并到一起
            output$score <- row_avg
            output$flag <- flag_tmp
            output$Ancestral_cell <- link[link$flag == flag_tmp,"Ancestral_cell"] # 添加假祖先cell type
            output$New_celltype <- link[link$flag == flag_tmp,"New_celltype"]  #添加new cell type
            output$Sister_celltype <- link[link$flag == flag_tmp,"Sister_celltype"]  #添加sister cell type
            sorted_df <- output[order(-output$score), ]
            new_df <- select(sorted_df,"flag","Ancestral_cell","New_celltype","Sister_celltype","GeneID","score","avg_log2FC","pct.1","pct.2","p_val_adj","S.avg_log2FC","S.pct.1","S.pct.2","S.p_val_adj")
            names(new_df)[names(new_df)=="GeneID"] <- "key_TF"
            results <- rbind(results,new_df)
        }
    }
}

write.table(results,file=paste0(output_dir,"/output/2.positive_NewCellTypes_key_TFs_score_stat.txt"),row.names = FALSE, sep="\t")


### ——step3: 预测celltype的负调控key TF和计算得分
### 负调控key TF
output2 <- data.frame()
for(i in 1:nrow(link)){
    if(is.na(link$Sister_celltype[i]) || link$Sister_celltype[i] == ""){
    	#message("sister_cell of the ", i, " row is null.") ## 有的new cell type没有sister cell type.
    	A <- read.table(paste0(output_dir,"/markers/N",i,"_vs_A",i,".cluster.AllmarkerGene.negative.txt"),sep="\t",header=T)
    	#S <- read.table(paste0("./markers/N",i,"_vs_S",i,".cluster.AllmarkerGene.txt"),sep="\t",header=T)
    	A_genes <- rownames(A)
    	#S_genes <- rownames(S)
    	#key_markers <- intersect(A_genes, S_genes)
    	len1 <- length(A_genes)
    	key_TF <- intersect(A_genes,TF)
    	keyTF_len <- length(key_TF)
    	print(keyTF_len)

    	prefix = paste0("N",i,"_vs_A",i)
    	key_TF <- paste(key_TF,collapse = ",")
    	result <- data.frame(Flag=prefix,Sister_celltype="No",Markers_num=len1,keyTF_num=keyTF_len,TF_IDs=key_TF)
    	output2 <- rbind(output2,result)

    }else{  
    	ancester <- read.table(paste0(output_dir,"/markers/N",i,"_vs_A",i,".cluster.AllmarkerGene.negative.txt"),sep="\t",header=T)
    	sister <- read.table(paste0(output_dir,"/markers/N",i,"_vs_S",i,".cluster.AllmarkerGene.negative.txt"),sep="\t",header=T)
    	#S <- read.table(paste0("./markers/N",i,"_vs_S",i,".cluster.AllmarkerGene.txt"),sep="\t",header=T)
    	ancester_genes <- rownames(ancester)
    	sister_genes <- rownames(sister)
    	key_markers <- intersect(ancester_genes, sister_genes)
    	print(length(ancester_genes))
    	print(length(sister_genes))
    	len1 <- length(key_markers)
    	print(len1)

    	key_TF <- intersect(key_markers,TF)
    	keyTF_len <- length(key_TF)
    	print(keyTF_len)

    	prefix = paste0("N",i,"_vs_A",i)
    	key_TF <- paste(key_TF,collapse = ",")
    	result <- data.frame(Flag=prefix,Sister_celltype="Yes",Markers_num=len1,keyTF_num=keyTF_len,TF_IDs=key_TF)
    	output2 <- rbind(output2,result)
    }
}
write.table(output2,file = paste0(output_dir,"/output/1.negative_NewCelltype_key_TFs_output.txt"), row.names = FALSE, sep="\t")

### 计算负调控key TF的score
results <- data.frame()
for (i in 1:nrow(output2)){
    if(output2$Sister_celltype[i] == "No"){  ## 有些new cell type是没有sister cell type的
        markers <- read.table(paste0(output_dir,"/markers/N",i,"_vs_A",i,".cluster.AllmarkerGene.negative.txt"),sep="\t",header=T)
        flag_tmp <- paste0("N",i,"_vs_A",i)
        markers$GeneID <- rownames(markers)
        TF_IDs <- output2[output2$Flag == flag_tmp, "TF_IDs"]
        TF_genes <- unlist(strsplit(TF_IDs, ","))
        #print(length(TF_genes))
        if(length(TF_genes)>0){  ## 有的new cell type没有key TFs.
            A_output <- markers[markers$GeneID %in% TF_genes, c("avg_log2FC","pct.1","pct.2","p_val_adj","GeneID")]
            temp_df <- as.data.frame(matrix(NA, nrow = nrow(A_output), ncol = ncol(A_output)))  # 生成一个空的data.frame，用于跟后面有sister>的情况合并；
            colnames(temp_df) <- paste("S", colnames(A_output), sep = ".")
            output <- cbind(A_output, temp_df)
            output$score <- scale(output$avg_log2FC)
            output$flag <- flag_tmp
            output$Ancestral_cell <- link[link$flag == flag_tmp,"Ancestral_cell"] # 添加假祖先cell type
            output$New_celltype <- link[link$flag == flag_tmp,"New_celltype"]  #添加new cell type
            #output$sister_cell <- link[link$flag == flag_tmp,"sister_cell"]  #添加sister cell type ——值为空
            output$Sister_celltype <- "No"
            sorted_df <- output[order(-output$score), ]
            #new_df <- select(sorted_df,"flag","ancestral_cell","New_celltype","GeneID","score","avg_log2FC","pct.1","pct.2","p_val_adj")
            new_df <- select(sorted_df,"flag","Ancestral_cell","New_celltype","Sister_celltype","GeneID","score","avg_log2FC","pct.1","pct.2","p_val_adj","S.avg_log2FC","S.pct.1","S.pct.2","S.p_val_adj")
            names(new_df)[names(new_df)=="GeneID"] <- "key_TF"
            results <- rbind(results,new_df)
        }
    }else{
        print(i)
        flag_tmp <- paste0("N",i,"_vs_A",i)
        ancester_markers <- read.table(paste0(output_dir,"/markers/N",i,"_vs_A",i,".cluster.AllmarkerGene.negative.txt"),sep="\t",header=T)
        sister_markers <- read.table(paste0(output_dir,"/markers/N",i,"_vs_S",i,".cluster.AllmarkerGene.negative.txt"),sep="\t",header=T)
        ancester_markers$GeneID <- rownames(ancester_markers)
        sister_markers$GeneID <- rownames(sister_markers)
        TF_IDs <- output2[output2$Flag == flag_tmp, "TF_IDs"]
        TF_genes <- unlist(strsplit(TF_IDs, ","))
        print(length(TF_genes))
        if(length(TF_genes)>0){  ## 有的new cell type没有key TFs.
            A_output <- ancester_markers[ancester_markers$GeneID %in% TF_genes, c("avg_log2FC","pct.1","pct.2","p_val_adj","GeneID")]
            S_output <- sister_markers[sister_markers$GeneID %in% TF_genes, c("avg_log2FC","pct.1","pct.2","p_val_adj","GeneID")]
            scaled_log_fold_change1 <- scale(A_output$avg_log2FC)  ## 对祖先avg_logFC进行scale
            scaled_log_fold_change2 <- scale(S_output$avg_log2FC)  ## 对sister细胞的avg_logFC进行scale
            matrix_data <- cbind(scaled_log_fold_change1, scaled_log_fold_change2)  # 合并2列;
            row_avg <- apply(matrix_data, 1, function(x) mean(x))  ## 求2个值的平均值。
            colnames(S_output) <- paste("S", colnames(S_output), sep = ".")  # 由于列名相同，在sister markers的列名上加上前缀“S.”
            output <- cbind(A_output, S_output)  ## 把祖先和sister的差异倍数等信息合并到一起
            output$score <- row_avg
            output$flag <- flag_tmp
            output$Ancestral_cell <- link[link$flag == flag_tmp,"Ancestral_cell"] # 添加假祖先cell type
            output$New_celltype <- link[link$flag == flag_tmp,"New_celltype"]  #添加new cell type
            output$Sister_celltype <- link[link$flag == flag_tmp,"Sister_celltype"]  #添加sister cell type
            sorted_df <- output[order(-output$score), ]
            new_df <- select(sorted_df,"flag","Ancestral_cell","New_celltype","Sister_celltype","GeneID","score","avg_log2FC","pct.1","pct.2","p_val_adj","S.avg_log2FC","S.pct.1","S.pct.2","S.p_val_adj")
            names(new_df)[names(new_df)=="GeneID"] <- "key_TF"
            results <- rbind(results,new_df)
        }
    }
}
write.table(results,file=paste0(output_dir,"/output/2.negative_NewCellTypes_key_TFs_score_stat.txt"),row.names = FALSE, sep="\t")


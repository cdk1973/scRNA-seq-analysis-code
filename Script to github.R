library(Seurat)
library(SeuratObject)
library(SeuratData)
library(stringr)
library(patchwork)
library(scales)
library(data.table)
library(reshape)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(ggsci)
library(ggpubr)
library(ggthemes)
library(ggrepel)
library(ggstream)
library(ggridges)
library(gghalves)
library(RColorBrewer)
library(cowplot)
library(pheatmap)
library(ggrastr)
library(future)
library(patchwork)
library(scales)
library(data.table)
library(reshape)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(ggsci)
library(ggpubr)
library(ggthemes)
library(ggrepel)
library(ggstream)
library(ggridges)
library(gghalves)
library(RColorBrewer)
library(cowplot)
library(pheatmap)
require(DOSE)
library(monocle)
library(ggpointdensity)
library(ggpmisc)
set.seed(5)
plan("multicore", workers = 25)
options(future.globals.maxSize = 300000 * 1024^2)#100000MB~=100G

#01 creat seurat object------------------------------------------------------------
##01.1 load data---------------------------------
for (i in seq(1:sample_number)) {
  #load data from the cellranger out
  data.data <- Read10X(paste0(mapping_dir1, name_list[i],"/outs/filtered_feature_bc_matrix"))
  tmp <- CreateSeuratObject(counts = data.data, min.cells = 3, min.features = 200, project = sample_list[i])
  tmp$sample <- name_list[i]
  tmp$disease<- dse_list[i]
  tmp[["percent.MT"]] <- PercentageFeatureSet(tmp, pattern = "^MT-")
  p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
  ggsave(paste0(qc_dir,name_list[i],"_mt_vln.pdf"), plot = p, width = 10, height = 6)
  ggsave(paste0(qc_total,name_list[i],"_mt_vln.png"), plot = p, width = 10, height = 6)
  
  plot1 <- FeatureScatter(tmp, feature1 = "nCount_RNA", feature2 = "percent.MT")
  plot2 <- FeatureScatter(tmp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- FeatureScatter(tmp, feature1 = "percent.MT", feature2 = "nFeature_RNA")
  plot4 <- plot1 + plot2+ plot3
  ggsave(paste0(qc_dir,name_list[i],"_mt_Scatter.pdf"),plot = plot4, width = 12, height = 4)
  ggsave(paste0(qc_total,name_list[i],"_mt_Scatter.png"),plot = plot4, width = 12, height = 4)
  
  #
  tmp[["percent.RBL"]] <- PercentageFeatureSet(tmp, pattern = "^RP[SL]")
  p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.RBL"), ncol = 3)
  ggsave(paste0(qc_dir,name_list[i],"_rbl_vln.pdf"), plot = p, width = 10, height = 6)
  ggsave(paste0(qc_total,name_list[i],"_rbl_vln.png"), plot = p, width = 10, height = 6)
  
  plot1 <- FeatureScatter(tmp, feature1 = "nCount_RNA", feature2 = "percent.RBL")
  plot2 <- FeatureScatter(tmp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- plot1 + plot2
  ggsave(paste0(qc_dir,name_list[i],"_rbl_Scatter.pdf"),plot = plot3, width = 8, height = 4)
  ggsave(paste0(qc_total,name_list[i],"_rbl_Scatter.png"),plot = plot3, width = 8, height = 4)
  
  #
  saveRDS(tmp, file = paste0(qc_dir,name_list[i],"_before_qc.rds"))
  assign(sample_list[i], tmp)
  print(paste0(sample_list[i],': ',ncol(tmp),' cells'))
}
#set mito_percent ########################
for (i in seq(1:sample_number)) {
  tmp <- get(sample_list[i])
  tmp <- subset(tmp, subset = percent.MT < 15 & nFeature_RNA > 200 & nFeature_RNA < 6000)
  p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
  ggsave(paste0(qc_dir,name_list[i],"_qc_vln.pdf"), plot = p, width = 10, height = 6,limitsize = F)
  ggsave(paste0(qc_total,name_list[i],"_qc_vln.png"),plot = p, width = 10, height = 6,limitsize = F)


  saveRDS(tmp, file = paste0(qc_dir,name_list[i],"_after_qc.rds"))
  tmp <- RenameCells(tmp, add.cell.id = name_list[i])
  assign(sample_list[i], tmp)
}
#
#normalization #####################
for (i in seq(1:sample_number)) {
  tmp <- get(sample_list[i])
  print(date())
  print(paste0(sample_list[i], ': SCTransform started'))
  tmp <- SCTransform(tmp, vars.to.regress = "percent.MT", verbose = FALSE,do.scale = T)
  print(date())
  print(paste0(sample_list[i], ': SCTransform finished'))
  #
  saveRDS(tmp, file = paste0(qc_dir,name_list[i],"_SCT.rds"))
  #
  tmp <- RunPCA(tmp, verbose=F)
  #
  png(paste0(qc_dir,name_list[i],"_pca_heatmap.png"), width=1000,height=2000)
  p=DimHeatmap(tmp, dims=1:30, cells=500, balanced=T)
  print(p)
  dev.off()
  #
  png(paste0(qc_dir,name_list[i],"_ElbowPlot.png"), height = 600, width = 700)
  p<- ElbowPlot(tmp, ndims = 30)
  print(p)
  dev.off()
  #
  saveRDS(tmp, file = paste0(qc_dir,name_list[i],"_bfPCR.rds"))
  print(paste0("--------------", sample_list[i], " completed (", i, "/" , sample_number ,")--------------"))
  ##PC_selection #####################
  tmp <- RunUMAP(tmp, dims = 1:20, verbose=F)
  tmp <- FindNeighbors(tmp, reduction = "pca", dims = 1:20)
  tmp <- FindClusters(tmp, res=res[i])
  tmp[["cluster"]] <- Idents(tmp)
  UMAP <- DimPlot(object = tmp, reduction = "umap", label = TRUE)
  ggsave(paste0(qc_dir,name_list[i],"_umap.pdf"), plot = UMAP, width = 8, height = 6,limitsize = F)
  ggsave(paste0(qc_dir,name_list[i],"_umap.png"), plot = UMAP, width = 8, height = 6,limitsize = F)
  saveRDS(tmp, file = paste0(qc_dir,name_list[i],"_PCR.rds"))
  assign(sample_list[i], tmp)
  print(paste0("--------------", sample_list[i], " completed (", i, "/" , sample_number ,")--------------"))
}
#Identify doublets ####
pK.df <- data.frame(matrix(nrow=0, ncol=2))
colnames(pK.df) <- c("Sample", "Optimal_pK")
#
for (i in c(1:sample_number)){
  sweep.res.list <- paramSweep_v3(get(sample_list[i]), PCs = 1:20, sct = T, num.cores=16)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- arrange(bcmvn, desc(BCmetric))$pK[1]
  tmp <- data.frame(Sample=name_list[i], Optimal_pK=pK)
  pK.df <- rbind(pK.df, tmp)
  print(bcmvn)
  print(paste0("--------------", sample_list[i], " completed (", i, "/" , sample_number ,")--------------"))
}
write.csv(pK.df,file =paste0(qc_dir,ts,"_Optimal_pK.csv"), sep = "" )
#
sample_inforation <- data.frame(matrix(nrow=0, ncol=6))
colnames(sample_inforation) <- c("Sample", "Number","Doublet_prop","AvailableCellNumber","Cell_num_pre","Ratio_use")
#
for (i in seq(1:sample_number)) {
  tmp <- get(sample_list[i])
  pK.use <- as.numeric(as.character(pK.df$Optimal_pK[i]))
  homotypic.prop <- modelHomotypic(tmp@meta.data$cluster)
  ##calculate ratio ##############################
  cell_num_pre <- as.numeric(length(tmp@meta.data$orig.ident))
  ratio <- cell_num_pre/1000*0.008
  ##################################################
  nExp_poi <- round(ratio*length(tmp@meta.data$orig.ident))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  tmp <- doubletFinder_v3(tmp, PCs = 1:20, pN = 0.25, pK = pK.use, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = T)
  tmp[["doublet"]] <- tmp[[paste("DF.classifications_0.25", pK.use, nExp_poi.adj, sep="_")]]
  prop <- nExp_poi.adj/length(tmp@meta.data$cluster)
  prop.tmp <- data.frame(Sample=name_list[i], Number=nExp_poi.adj, Doublet_prop=prop,
                         AvailableCellNumber=length(tmp@meta.data$doublet[tmp@meta.data$doublet=='Singlet']),
                         Cell_num_pre=cell_num_pre,Ratio_use=ratio)
  #
  sample_inforation <- rbind(sample_inforation, prop.tmp)
  saveRDS(tmp, file = paste0(qc_dir,name_list[i],"_doublets.rds"))
  assign(sample_list[i], tmp)
  print(paste0("--------------", sample_list[i], " completed (", i, "/" , sample_number ,")--------------"))
}
#
write.csv(sample_inforation,file =paste0(qc_dir,ts,"_sample_info.csv"), sep = "" )
#
for (i in seq(1:sample_number)){
  doub <- DimPlot(get(sample_list[i]), group.by='doublet', cols=c('firebrick', 'grey90'))
  pdf(paste0(qc_dir,name_list[i],"_UMAP_doublet.pdf"), height = 6, width = 8)
  print(doub)
  dev.off()
  png(paste0(qc_dir,name_list[i],"_UMAP_doublet.png"), height = 600, width = 800)
  print(doub)
  dev.off()
}
#filter out doublet ###############
for (i in seq(1:sample_number)){
  tmp <- get(sample_list[i])
  tmp <- subset(tmp, doublet=='Singlet')
  tmp <- SCTransform(tmp, verbose = FALSE, do.scale = T)
  saveRDS(tmp, file = paste0(qc_dir,name_list[i],"_final.rds"))
  print(paste0(i, ' completed', " (", i, "/" , sample_number ,")"))
}

#Integration ###########################
for (i in c(1:sample_number)) {
  tmp <- readRDS(file = paste0(qc_dir,sample_list[i],"_final.rds"))
  assign(sample_list[i],tmp)
}

int.list <- c()
for (i in sample_list) {
  tmp <- get(i)
  tmp2 <- list(tmp)
  int.list <- c(int.list,tmp2)
}

saveRDS(int.list,paste0(analysis_dir,'int.list.rds'))

int.list <- readRDS(paste0(analysis_dir,'int.list.rds'))
int.features <- SelectIntegrationFeatures(object.list = int.list, nfeatures = 3000)
int.list <- PrepSCTIntegration(object.list = int.list, anchor.features = int.features, verbose = FALSE)

int.anchors <- FindIntegrationAnchors(object.list = int.list, normalization.method = "SCT",
                                      anchor.features = int.features, verbose = FALSE)
saveRDS(int.anchors, paste0(analysis_dir,ts,'_int_anchors2.rds'))


combination <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT", verbose = FALSE)
DefaultAssay(combination) <- "integrated"
saveRDS(combination, paste0(analysis_dir,ts,'_int_seu2.rds'))

#02 Cell identification--------------
pc.list <- pc.num <- c(1:22)
n.neighbors=10
min.dist=0.3
spread=1.1
resolution=1
combination <- FindNeighbors(combination, reduction = "pca",dims = pc.list,  assay= "integrated")
combination <- FindClusters(object = combination, resolution = resolution)
combination <- RunUMAP(combination, reduction = "pca",dims = pc.list,  assay= "integrated",
                       n.neighbors= n.neighbors,
                       min.dist= min.dist,
                       spread= spread
)

celltypes <- c('Granulosa cells (GC)',
               'Stromal cells (SC)',
               'Ovarian surface epithelial cells (OSEC)',
               'Theca cells (Theca)',
               'Endothelial cells (EC)',
               'Fibroblasts (Fib)',
               'Macphages type 1 (Mac1)',
               'Macphages type 2 (Mac2)',
               'Neutrophils (Neu)'
               
)


C1 <- c(8,1,4,15,3,2,6,0,19)
C2 <- c(7,13,23,10)
C3 <- c(11)
C4 <- c(20)
C5 <- c(5,16,14)
C6 <- c(9,18,21)
C7 <- c(12)
C8 <- c(22)
C9 <- c(17)


df_cl_clt <- NULL
for (i in 1:length(celltypes)) {
  m <- paste0('C',i)
  tmp <- get(m)
  df_tmp <- data.frame(cell=celltypes[i],cluster=tmp)
  df_cl_clt <- rbind(df_cl_clt,df_tmp)
}
new.cluster.ids <- df_cl_clt$cell
names(new.cluster.ids) <-  df_cl_clt$cluster

combination <- RenameIdents(combination, new.cluster.ids)
Idents(combination) <- factor(Idents(combination), levels=celltypes)

combination$celltype <- Idents(combination)
UMAP <- DimPlot(combination, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 6,
                cols= brewer.pal(13,"Paired")[1:9]
)
ggsave(paste0(clusters_dir,ts,"_celltype.png"), plot = UMAP, width = 12, height = 8)
ggsave(paste0(clusters_dir,ts,"_celltype.pdf"), plot = UMAP, width = 12, height = 8)



combination2 <- combination
#combination <- combination2

cell_abbr <- c('GC',
               'SC',
               'OSEC',
               'Theca',
               'EC',
               'Fib',
               'Mac1',
               'Mac2',
               'Neu'
)
names(cell_abbr) <- levels(combination)
combination <- RenameIdents(combination, cell_abbr)

Idents(combination) <- factor(Idents(combination), levels=cell_abbr)
combination$cell_abbr <- Idents(combination)

UMAP <- DimPlot(combination, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 6,
                cols= brewer.pal(13,"Paired")[1:9]
)
ggsave(paste0(clusters_dir,ts,"_cell_abbr.pdf"), plot = UMAP, width = 10, height = 8)
ggsave(paste0(clusters_dir,ts,"_cell_abbr.png"), plot = UMAP, width = 10, height = 8)

UMAP <- DimPlot(combination, reduction = "umap", label = F, pt.size = 0.5,label.size = 6,
                cols= brewer.pal(13,"Paired")[1:9])
ggsave(paste0(clusters_dir,ts,"_cell_abbr_noLab.png"), plot = UMAP, width = 10, height = 8)
saveRDS(combination,paste0(ana_dir,ts,"_celltype.rds"))

#03 cellmarker----------------------------------------------------
DefaultAssay(combination) <- "RNA"
df <- FindAllMarkers(combination)
df$gene <- rownames(df)
df <- df[df$avg_log2FC>0.5&df$p_val_adj<=0.05,]
df <- df%>%arrange(cluster,desc(avg_log2FC))
write.csv(df,paste0(clusters_dir,'cellmrk.csv'))

df <- read.csv(paste0(clusters_dir,'cellmrk.csv'))
df$cluster <- factor(df$cluster,levels = c('GC',
                                           'SC',
                                           'OSEC',
                                           'Theca',
                                           'EC',
                                           'Fib',
                                           'Mac1',
                                           'Mac2',
                                           'Neu'))
df2 <- df[df$avg_log2FC>0.25&df$p_val_adj<=0.05,]
write.csv(df2,paste0(clusters_dir,'marker_gn_tt.csv' ),row.names = F)
cellmarkers <- df2 
cellmarkers$gene <- apply(cellmarkers,1,function(vec){strsplit(vec['gene'],'\\.')[[1]][1] })

mark <- cellmarkers[cellmarkers$p_val_adj<0.01&cellmarkers$avg_log2FC>0.25,]
mark <- mark[order(mark$avg_log2FC,decreasing = T),]
gene_ls <- c()
for (cell in levels(mark$cluster)) {
  tmp <- mark[mark$cluster==cell,]
  gene_tmp <- tmp$gene[1:50]
  gene_ls <- c(gene_ls,gene_tmp)
}
write.csv(gene_ls,paste0(clusters_dir,"top_mkr.csv"),row.names = F)
DefaultAssay(combination) <- "RNA"
cluster.averages <- AverageExpression(combination, return.seurat = TRUE)
mtx <- cluster.averages@assays$RNA@scale.data[gene_ls,]
hmcols <- c(colorRampPalette(c('#8B13A8', '#F4F4F4'))(10),
            colorRampPalette(c('#F4F4F4', "#ff4800"))(10)
)
pheatmap::pheatmap(mtx,
                   cluster_rows = F,
                   cluster_cols =F,
                   scale = "row",#"none", "row", "column"
                   color = hmcols,
                   show_rownames = F,
                   border_color = NA, #"grey60",
                   fontsize = 14,
                   fontsize_row = 10, #fontsize,
                   fontsize_col = 20, #fontsize,
                   fontsize_number = 0.8* fontsize,
                   kmeans_k = NA,
                   cutree_rows = NA,
                   treeheight_row=0,
                   cutree_cols = NA,
                   legend_breaks = NA,
                   legend_labels = NA,
                   annotation = NA,
                   annotation_col= NA,
                   annotation_legend = TRUE,
                   drop_levels = TRUE,
                   #gaps_row = seq(1,nrow(heatmap_up)),
                   gaps_col = c(1:9),
                   gaps_row = c(seq(0,450,50)),
                   labels_row = NULL,
                   labels_col = NULL,
                   cellwidth = 30,
                   cellheight = 2,
                   revC=F,
                   # # width = wd,
                   # # height = ht,
                   # fontface="italic",
                   # main='',
                   filename=paste0(clusters_dir,'mrk_htmp.pdf')
)



#04 DEG analysis after deconvolution--------------------------
S2_S1_tmp <- NULL
for (cell in unique(seurat_rds@meta.data$cell_name) ) {
  tmp <- FindMarkers(seurat_rds, ident.1=paste0(cell, '_Disease'), ident.2=paste0(cell, '_Healthy'),logfc.threshold = 0.2)
  tmp$gene <- rownames(tmp)
  tmp$celltype <- cell
  tmp$group <- paste0('Disease_Healthy')
  S2_S1_tmp <- rbind(S2_S1_tmp, tmp)
  print(paste0(cell, ' is finished'))
}
S2_S1_tmp <- S2_S1_tmp[S2_S1_tmp$p_val_adj<0.05&abs(S2_S1_tmp$avg_log2FC)>0.25,]
write.csv(S2_S1_tmp,paste0(DEG_dir,'Dse_Hty.csv'),row.names = F)

df <- read.csv(paste0(DEG_dir,'Dse_Hty.csv'))
df$gene <-  toupper(df$gene)

df$drc <- ifelse(df$avg_log2FC>0,'U','D')
for (pat in c('U','D')) {
  df_pat <- df[df$drc==pat,]
  df_GO <- NULL
  for (cell in unique(df_pat$celltype)) {
    tmp <- df_pat[df_pat$celltype==cell,]
    gene_tmp <- tmp$gene
    tmp <- data.frame(gn=gene_tmp)
    colnames(tmp) <- cell
    tmp_gene <- data.frame(t(tmp))
    df_GO <- plyr::rbind.fill(df_GO,tmp_gene)
  }
  df_GO <- t(df_GO)
  colnames(df_GO) <- unique(df_pat$celltype)
  df_GO[is.na(df_GO)] <- ''
  write.csv(df_GO,paste0(DEG_dir,pat,'_forGO.csv'),row.names = F)
}


for (drc in c('U','D')) {
  df_GO <- read.csv(paste0(DEG_dir,'GO_',drc,'.csv'))
  df_GO <- melt(df_GO,id.vars = 'Description',variable.name ='celltype',value.name='logP' )
  if (drc=='U') {
    cor_slt <- '#da2d20'
  }else {cor_slt <- '#33b2e8'}
  dat2 <- df_GO
  #dat2$logP <- -dat2$logP
  dat2$celltype <- factor(dat2$celltype,levels = c('GC',
                                                   'SC',
                                                   'OSEC',
                                                   'Theca',
                                                   'EC',
                                                   'Fib',
                                                   'Mac1',
                                                   'Mac2',
                                                   'Neu') )
  dat3 <- cast(dat2,Description~celltype,value = 'logP')
  dat3[is.na(dat3)] <- 0
  dat4 <- as.data.frame(dat3)
  rownames(dat4) <- dat4$Description
  dat4[dat4>0] <- 1
  dat4[dat4<0] <- -1
  dat4 <- dat4%>%dplyr::mutate(pat_num=rowSums(as.data.frame(lapply(dplyr::select(dat4,2:c(ncol(dat4)) ), as.numeric)))   )
  identical(dat3$Description,rownames(dat4))
  dat3$pat_num1 <- dat4$pat_num
  dat3 <- dat3%>%dplyr::mutate(pat_num=rowSums(dplyr::select(.,2:c(ncol(dat3)-1) )),
                               sort_order = purrr::pmap_dbl(dplyr::select(.,2:c(ncol(dat3)-1)   ), ~ { sum(2^(length(c(...)) - which(c(...) != 0))) })
  )%>%
    arrange( desc(pat_num1),-sort_order,desc(pat_num))
  dat3 <- as.data.frame(dat3)
  dat2$Description <- factor(dat2$Description,levels = rev(dat3$Description))
  ggplot(dat2,aes(x=celltype,y=Description,color=logP))+
    geom_point(size=5)+
    scale_color_gradient2(low ='white' ,
                          high = cor_slt)+
    theme(panel.border = element_rect(fill = 'transparent',colour ='black',size=1),
          panel.grid.minor = element_blank(),
          panel.grid.major=element_blank(),
          panel.background = element_rect(fill = NA,colour = NA),
          # plot.background = element_rect(fill = 'transparent',colour = NA),
          # axis.line.x = element_line(size=0.2, colour = "black"),
          axis.ticks= element_line( colour = "black"),
          # axis.ticks.x = element_line(size=0.2, colour = "black"),
          #axis.text.y = element_blank(),
          axis.text.x=element_text(),
          # axis.title.x = element_text(size=8,color="black",face="plain"),
          axis.title=element_blank(),
          # #legend.text=element_blank(),
          # legend.title=element_blank(),
          #legend.position = "none",
          # plot.title = element_text(hjust = 0.5,size=12),4
          
          plot.margin=unit(c(0,0,0,0),'lines'))
  
  ggsave(paste0(DEG_dir,'GO_',drc,'.pdf'),width = 6,height = 4.5)
  
}



#05 monocle analysis--------------------------------------------
mncl <- subset(seurat_rdscell_name=='SC')
exp_mat <- as.matrix(GetAssayData(mncl, slot='data'))
exp_mat <- exp_mat[rowSums(exp_mat)!=0,]
pd <- data.frame(mncl[[]])
fd <- data.frame(gene_short_name=rownames(exp_mat))
rownames(fd) <- fd$gene
pd <- new("AnnotatedDataFrame", data = pd)
fd <- new("AnnotatedDataFrame", data = fd)

mncl= FindVariableFeatures(mncl)
varGenes = VariableFeatures(mncl)

ordering_genes = head(varGenes, 1000)
cds_obj <- newCellDataSet(as(exp_mat, "sparseMatrix"),
                          phenoData = pd,
                          featureData = fd)
cds_obj <- estimateSizeFactors(cds_obj)
cds_obj <- estimateDispersions(cds_obj)
cds_obj <- setOrderingFilter(cds_obj, ordering_genes)
cds_obj <- reduceDimension(cds_obj, max_components = 2,
                           method = 'DDRTree')
cds_obj <- orderCells(cds_obj)

saveRDS(cds_obj,paste0(monocle_dir,"Epithelial_monocle.rds"))


p=plot_cell_trajectory(cds_obj,color_by="Pseudotime")+
  scale_color_gradientn(colours = scales::viridis_pal(option = "D")(20))
ggsave(paste0(monocle_dir,'psedutime_traj.png'),plot=p,width = 4,height = 3.5)
ggsave(paste0(monocle_dir,'psedutime_traj.pdf'),plot=p,width = 4,height = 3.5)

p=plot_cell_trajectory(cds_obj,color_by = "spp")+
  facet_wrap(spp~.,ncol=1)+
  scale_color_manual(values= c("#2C77E5", "#E0760B") )
ggsave(paste0(monocle_dir,'psedutime_dese.png'),plot=p,width = 4,height = 6)
ggsave(paste0(monocle_dir,'psedutime_dese.pdf'),plot=p,width = 4,height = 6)


plotdf=data.frame(Pseudotime=pData(cds_obj)$Pseudotime,cluster=pData(cds_obj)$cell_new,state=pData(cds_obj)$State,age=pData(cds_obj)$spp)
plotdf=subset(plotdf,plotdf$Pseudotime != 'Inf')
write.csv(plotdf,paste0(monocle_dir,'pse_time_prop.csv'),row.names = F)

plotdf <- read.csv(paste0(monocle_dir,'pse_time_prop.csv'))
plotdf$age <- factor(plotdf$age,levels = c('Normal','Tumor'))
plotdf$state <- factor(plotdf$state,levels = c('1','2','3'))
ggplot(plotdf, aes(x=Pseudotime,y=state,fill=age))+
  geom_density_ridges(alpha=1,#transparent degree
                      scale=2,panel_scaling=F,
                      #color='black',
                      na.rm=TRUE,
                      quantile_lines=TRUE,
                      vline_linetype=2,
                      vline_color='black',
                      quantile_fun=function(x,...)mean(x)
                      
  ) +
  scale_y_discrete("")+
  scale_fill_manual(values=c("#2C77E5", "#E0760B"))+
  theme_minimal()+ 
  geom_rug(aes(color=plotdf$age))+
  theme(panel.grid = element_blank())
ggsave(paste0(monocle_dir,'psedutime_prop.png'),plot=p,width = 4,height = 3)
ggsave(paste0(monocle_dir,'psedutime_prop.pdf'),plot=p,width = 4,height = 3)

BEAM_res=BEAM(cds_obj,branch_point = 1,cores = 20)
BEAM_res=BEAM_res[,c("gene_short_name","pval","qval")]

tmp1=plot_genes_branched_heatmap(cds_obj[row.names(subset(BEAM_res,qval<1e-4)),],
                                 branch_point = 1,
                                 num_clusters = 4, 
                                 cores = 20,
                                 branch_labels = c("Cell fate 1", "Cell fate 2"),
                                 branch_colors = c("#979797", "#2C77E5", "#E0760B"), 
                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 return_heatmap = T 
)

pdf(paste0(monocle_dir,"branched_heatmap.pdf"),width = 5,height = 6)
p <- tmp1$ph_res
print(p)
dev.off()






































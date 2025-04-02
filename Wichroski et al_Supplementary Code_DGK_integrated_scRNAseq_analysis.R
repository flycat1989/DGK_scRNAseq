rm(list=ls())
library(Seurat)
library(ggpubr)


load('/integrated_obj_sele.RData')

pdf('/umap_CLUSTER.pdf',width = 10,height = 8)
print(DimPlot(obj_sele, reduction = "umap",label = TRUE))
dev.off()
LabelClusters(plot = plot, id = "ident")

FeaturePlot(obj_sele,features = c('DGKZ','DGKA'),min.cutoff = 1, max.cutoff = 3,order = TRUE)

for(g in c('DGKA','DGKZ','PDCD1','TCF7','CD8A','TOX')){
  pdf(paste0('/umap_genes_',g,'.pdf'))
  print(FeaturePlot(obj_sele,features = g,cols = c("white", "#BE2BBB"),
                    max.cutoff = 3,sort.cell = TRUE,raster=FALSE
  ))
  dev.off()}


DoHeatmap(res1, features = c('DGKZ','DGKA'), size = 4,
          angle = 90) #+ NoLegend()

Idents(obj_sele) = 'Cell_Cluster_level2'

cts_sele = c("CD8+ Tem","NK",'Naive T')
res = subset(obj_sele, idents = cts_sele, invert = FALSE)

obj_sele@meta.data$Cell_Cluster_level2 = as.factor(obj_sele@meta.data$Cell_Cluster_level2 )
VlnPlot(res, features = c('DGKZ','DGKA','PDCD1'),cols=c('red','blue','purple'))#, split.by = "groups")

RidgePlot(res, features = c('DGKZ','DGKA','PDCD1'), ncol = 2)

res$Cell_Cluster_level2
VlnPlot(res, features = c('DGKZ','DGKA','PDCD1'))+ scale_colour_manual(values = c("red", "blue", "purple"))

# IFNG, GZMB, GZMA, GZMK, PDCD1, CTLA4, HAVCR2 (TIM3), LAG3, TIGIT, TOX, TCF7

DotPlot(res, features = c('DGKZ','DGKA','PDCD1')) + RotatedAxis()

#DoHeatmap(subset(pbmc3k.final, downsample = 100), features = features, size = 3)


expl = c()
ctl = c()
g_l = c()
for(i in c('DGKZ','DGKA','PDCD1')){
  expl = c(expl,res@assays$RNA@data[i,])
  ctl = c(ctl,res@meta.data[colnames(res@assays$RNA@data),'Cell_Cluster_level2'])
  g_l = c(g_l,rep(i,times=length(res@assays$RNA@data[i,])))
}
data_plot = data.frame("Expression level"=expl,
                       'Cell type'=ctl,
                       'gene'=g_l)
#data_plot = as.data.frame(t(as.data.frame(res@assays$RNA@data[c('DGKZ','DGKA','PDCD1'),])))
#data_plot[,'cell type'] = res@meta.data[row.names(data_plot),'Cell_Cluster_level2']
pdf('/markers_celltype_plot.pdf')
print(ggboxplot(data_plot, x = 'gene', y = "Expression.level",
                color = "Cell.type", palette =c('coral1','cornflowerblue','purple'),
                add = "jitter"))
dev.off()

pdf('/markers_celltype_plot_violin.pdf')
print(ggviolin(data_plot, x = 'gene', y = "Expression.level",
               color = "Cell.type", palette =c('coral1','cornflowerblue','purple'),
               add = "jitter"))
dev.off()



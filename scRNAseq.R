##### scRNA-seq assignment for gene expr analysis course #####

setwd("~/Documents/R for gene expression analysis/scRNAseq")
install.packages("dplyr")
install.packages("Seurat")
library(dplyr)
library(Seurat)
library(ggplot2)

#Read data
data <- Read10X("./mm10") #Read all relevant files to the current directory
#27998 genes and 9128 cells


##### Filtering of data #####

#Create a Seurat object where counts and all metadata will be stored using createSeuratObject(). 
#When creating the object, remove genes that are expressed in fewer than 3 cells, and remove cells that have 
#less then 200 genes expressed. 
so <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)
#17311 genes and 9084 cells, so 10687 genes were removed as well as 44 cells

#Calculate the percentage of mitochondrial gene expression for every cell by using PercentageFeatureSet(). 
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-") #Adds an extra element to the "so" structure

#Plot the 
#1) number of genes (nFeature_RNA), 
#2) number of counts (nCount_RNA) per cell and 
#3) mt-percentage (percent.mt) per cell using VlnPlot(). 
#If you want to plot all in the same plot then specify the column names as a list 
#features = c('colname','colname','colname') and add ncol = 3 in the end.

VlnPlot(so, features = c('nFeature_RNA','nCount_RNA','percent.mt'))

#Now we will try to decide which cells that might be duplicates or destroyed and needs to be filtered out. 
#Make two new plots:

FeatureScatter(so, "nCount_RNA", "percent.mt") #percent.mt <10 (nCount > 5000), nCount <40000
FeatureScatter(so, "nCount_RNA", "nFeature_RNA") #Use nFeat >800, nCount < 40000

so_clean <- so
so_clean <- subset(so_clean, subset = nCount_RNA < 25000)
so_clean <- subset(so_clean, subset = nFeature_RNA > 800)
so_clean <- subset(so_clean, subset = percent.mt < 10)
#Same number of genes left but only 8894 cells, 190 cells were removed

FeatureScatter(so_clean, "nCount_RNA", "percent.mt")
FeatureScatter(so_clean, "nCount_RNA", "nFeature_RNA")


##### Normalization and scaling #####

so_norm <- NormalizeData(so_clean) #Log-normalization is default
#Select genes with most cell-to-cell variation:
sel <- FindVariableFeatures(so_norm) #vst selection method is default
top10 <- head(VariableFeatures(sel), 10) #list of 10 top genes
VariableFeaturePlot(sel)
plot1 <- VariableFeaturePlot(sel)
LabelPoints(plot = plot1, points = top10)

sel_scaled <- ScaleData(sel) #Scaling before PCA

#PCA reduces number of dimensions:
sel_pca <- RunPCA(sel_scaled) # @reductions added to the structure, with "pca" under it
print(sel_pca[["pca"]], dims = 1:8, nfeatures = 2) #variables from excercise
#Dimension plot
DimPlot(sel_pca) #How should this be analyzed?
#Heatmaps for the first 9 pcas, use 500 cells
DimHeatmap(sel_pca, dims = 1:9, cells = 500)
#How many PCs should be used? Make Elbow plot:
ElbowPlot(sel_pca)
#I would probably use 14 PCs. If less, then 8.

##### Clustering cells #####

sel_snn <- FindNeighbors(sel_pca, dims = 1:14)
sel_clust_06 <- FindClusters(sel_snn, resolution = 0.6)  
sel_clust_08 <- FindClusters(sel_snn, resolution = 0.8) 
sel_clust_10 <- FindClusters(sel_snn, resolution = 1.0) 
#With tsne, the resolution needs to be 0.8 to separate clusters well (tried many)
#With umap, the resolution 0.6 is ok, with 14 PCs:
sel_umap <- RunUMAP(sel_clust_06, dims = 1:14)
sel_tsne <- RunTSNE(sel_clust_08, dims = 1:14)
plot_umap <- DimPlot(sel_umap, reduction = "umap")
plot_tsne <- DimPlot(sel_tsne, reduction = "tsne")

##### Differential expression analysis #####

#Find markers that separate two groups:
markers1_2 <- FindMarkers(sel_umap, ident.1 = 1, ident.2 = 2)
head(markers1_2)

#Find markers that separate two groups with several clusters in each:
markers12_34 <- FindMarkers(sel_umap, ident.1 = c(1,2), ident.2 = c(3,4))
head(markers12_34)

#Find markers that are specific for each cluster, compared to all other clusters:
neuron.markers <- FindAllMarkers(sel_umap, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4)
#Two/ten top markers for every structure:
markers_top2 <- neuron.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
markers_top2
markers_top10 <- neuron.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

#Visualize the expression of some of these genes using VlnPlot(), FeaturePlot() and DotPlot().
VlnPlot(sel_umap, features = markers_top2$gene[1]) #Testing function
for (i in c(1:4,13:18,25:28)) {
  vln_plot <- VlnPlot(sel_umap, features = markers_top2$gene[i]) 
  plot(vln_plot)
  assign(paste0("vln_plot", i), vln_plot)
} #Generates a vln plot for some of the top 2 genes of the 15 cell clusters

for (i in c(1:4,13:18,25:28)) {
  feat_plot <- FeaturePlot(sel_umap, features = markers_top2$gene[i]) 
  plot(feat_plot)
  assign(paste0("feat_plot", i), feat_plot)
} #Generates a feat plot for the same genes as the vln plots

for (i in c(1:4,13:18,25:28)) {
  dot_plot <- DotPlot(sel_umap, features = markers_top2$gene[i]) 
  plot(dot_plot)
  assign(paste0("dot_plot", i), dot_plot)
} #Generates a dot plot for the same genes as the vln plots

#Make a heatmap with top 10 genes for each cluster
heatmap_top10 <- DoHeatmap(sel_umap, features = markers_top10$gene)
#Warning message:
#  In DoHeatmap(sel_umap, features = markers_top10$gene) :
#  The following features were omitted as they were not found in the scale.data slot for the RNA assay: Ppp2r2b




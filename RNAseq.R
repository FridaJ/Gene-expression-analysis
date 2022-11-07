##### RNA-seq assignment for gene expr analysis course #####

setwd("Documents/R for gene expression analysis/RNA-seq/")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install() #To install latest version if multiple versions are installed
BiocManager::install(c("DESeq2", "pheatmap", "RColorBrewer","biomaRt"))
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(biomaRt)

#Read data
counts <- read.table("geneCounts.txt", header = TRUE, row.names = 1, as.is = TRUE)
info <- read.table("sample_info.txt", header = TRUE, row.names = 1, as.is = TRUE)
condition <- info$Group

dds <- DESeqDataSetFromMatrix(counts, data.frame(condition), ~condition) #A structure to be filled with processing info
assay(dds) #Shows the data in dds


##### Quality control of data before analysis #####

#There are two main data transformation functions within DESeq, one is the rlog (which tranforms the count data 
#to the log2 scale) and the other is vst (which calculates the variance stabilizing transformation from the 
#fitted dispersion-mean relation).

### Q:  When to use which? See report. rlog for machine learning!

rld <- rlogTransformation(dds)
vsd <- varianceStabilizingTransformation(dds)

#PCA plots:
pca_rld <- plotPCA(rld)
pca_vsd <- plotPCA(vsd) #Very similar to pca_rld

#Answers to report questions:
#It is easy to see the difference between the ctrl and VL groups. 
#67% of the data variance can be explained by the first two principles components.
#The number of genes used for the pca plots are 500 as set by default.

###Distance plots using pheatmap:

#Transpose data (use assay on transformed data)
t_rld <- t(assay(rld))
d <- dist(t_rld)
dm <- as.matrix(d)

#Make heatmap
annotation_col<-data.frame(Group=condition)
row.names(annotation_col) <- row.names(info) #Set row names of annotation_col to the sample names
annotation_colors <- list(Group = c(Ctrl="cadetblue2",VL="deeppink3"))

pheatmap(dm, cluster_rows = TRUE, show_rownames = TRUE, main = "Heatmap for log transformed count data")
pheatmap(dm, cluster_rows = TRUE, annotation_col = annotation_col, annotation_colors = annotation_colors, 
         main = "Heatmap for log transformed count data")

#If you had more metadata that you would like to add as annotation to your graph, you can just add the annotation 
#and the corresponding colors using these commands. Let's assume you have the sex of the samples:

#Sample  Sex
#Ctrl1   F
#Ctrl2   M
#Ctrl3   F
#Ctrl4   M
#Ctrl5   F
#VL1     M
#VL2     F
#VL3     M


##### Data analysis #####

#Use DESeq to analyze data:
dds <- DESeq(dds)
res <- results(dds)

#Histogram on p.adj values from results(dds) to view enrichment
hist(res$padj, main = "Histogram of adjusted p-values vs frequency", xlab = "Adjusted p-value",
     breaks = 20, col = "sky blue")
hist(res$padj, main = "Histogram of adjusted p-values vs frequency", xlab = "Adjusted p-value",
     breaks = 100, col = "violet")
#It can be seen that there are almost 2500 genes that have an adjusted p-value ≤ 0.05, whereof over 1200
#have an adjusted p-value ≤ 0.01

#MA-plot on data
plotMA(dds, ylim = c(-8, 8)) #Default ylim was 4-5 ish. Triangles means that points are outside the graph.
#THe highest and lowest fold changes were around 7.5 and -7.5, respectively, according to the plot.

#Remove genes with very low counts (p.adj = na). These can't be used for DE.
na_i <- which(is.na(res$padj) == TRUE) #Length is 6660, remove these rows from data:
res_no_na <- results(dds[-na_i,])
nrow(res_no_na) #19511 genes with high enough counts for DE


##### Getting and presenting results #####

#Which, and how many, genes have an adjusted p-value ≤ 0.01?
sign_i <- which(res_no_na$padj <= 0.01) #1227 significant genes
res_sig <- res_no_na[sign_i,]

#Sort genes on fold change:
gene_order <- order(res_sig$log2FoldChange, decreasing = TRUE) #Indices of genes in fold change decreasing order
res_sig_sorted <- res_sig[gene_order,] #Results with all significant genes in order of decreasing fold change
#The most up-regulated gene is ENSG00000174358 and the most down-regulated gene is ENSG00000142408.

#Sort genes on adjusted p-values:
gene_order_padj <- order(res_sig$padj) #Indices of genes in decreasing order of p adj
res_sig_sorted_padj <- res_sig[gene_order_padj,] #Results with all significant genes in order of decreasing padj
#The most significant gene (ENSG00000116133) is up-regulated (log2FC = 2.6)

#At this point, you can save your results into a file. You could use the write.table() function for instance.
write.table(res_sig_sorted_padj, file = "RNA_seq_results_padj.txt", sep = "\t", eol = "\r") #doesn't work?

### Plotting heatmaps ###

sel <- rownames(assay(rld))%in%rownames(res_sig) #Indices of sign genes in assay(rld) connected by row name in res_sig
sig2heatmap<-assay(rld)[sel,] #Genes chosen for heatmap
sig2heatmap_diff <- sig2heatmap - rowMeans(sig2heatmap) #Calculate the difference to the mean of each row
pheatmap(sig2heatmap_diff, show_rownames = FALSE, annotation_col = annotation_col, 
         annotation_colors = annotation_colors, main = "Heatmap for significant genes vs samples")

#From the heatmap it can be seen that the two groups VL and Ctrl each have a set of mainly up- vs down-regulated genes. 
#The genes are significant because when they are up-regulated for one group, they are down-regulated for the other, and vice versa.

#plotCounts for specific genes:
plotCounts(dds, "ENSG00000174358")
plotCounts(dds, "ENSG00000142408")

##### Optional tasks #####

###Select rlog transformed genes that correspond to 20 most significant genes (p adj)
genes_20 <- rownames(assay(rld))%in%rownames(res_sig_sorted_padj[1:20,]) #Indices of sign genes in assay(rld) connected by row name in res_sig
genes2heatmap<-assay(rld)[genes_20,] #Genes chosen for heatmap
genes2heatmap_diff <- genes2heatmap - rowMeans(genes2heatmap) #Calculate the difference to the mean of each row
pheatmap(genes2heatmap_diff, show_rownames = TRUE, annotation_col = annotation_col, 
         annotation_colors = annotation_colors, main = "Heatmap for 20 most significant genes")

### Annotation

#First specify which BioMart database we will be using with the useMArt() function. In our case select ensembl.
mart <- useMart("ensembl")
#Then set the dataset with the useDataset() function. Specify the dataset to hsapiens_gene_ensembl and the mart to ensembl.
mart <- useDataset("hsapiens_gene_ensembl", mart = mart)
#Now let's map our Ensembl IDs to their corresponding gene symbol using the getBM() function. 
#Use the attributes, filters, values and mart arguments. Store the results in an object.
anno_genes <- row.names(genes2heatmap)
anno <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"), 
              values = anno_genes, mart = mart)

#ERROR MESSAGE when using attribut GO_id:
#Error in curl::curl_fetch_memory(url, handle = handle) : 
#  transfer closed with outstanding read data remaining

#Finally, add this information to the results object you used in the previous step.









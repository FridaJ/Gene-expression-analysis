##### Microarray excercises for gene expr analysis course #####

setwd("Documents/R for gene expression analysis/Microarrays")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("estrogen", "limma", "affy","oligo","annotate","hgu95av2.db"))
library("estrogen") 
library("limma")
library("affy")
library("oligo")
library("annotate")
library("hgu95av2.db")

install.packages("gplots") # Because it said so in the task
library(gplots)
install.packages("ggplot2") # Because I like it :)
library(ggplot2)

readdir <- system.file("extdata", package = "estrogen") #Locate folder "extdata"
dir(readdir) #To view contents of library
targets<-readTargets("targLimma.txt", path=readdir)

### Reading the Affymetrics data ###

install.packages("stringr") #For regex
library(stringr)

rawData <- read.celfiles(filenames = file.path(readdir,targets$FileName))

#"By default the samples get the same names as the data files, but since you have 
#a sample name vector in the targets data frame you can use this together with 
#the sampleNames argument to give the samples proper names."

raw_data_rma <- rma(rawData, background=FALSE, normalize=FALSE) #Default is bgcorr=TRUE!
raw_data_rma_bg <- rma(rawData, background=TRUE, normalize=FALSE)
raw_data_rma_norm <- rma(rawData, background=FALSE, normalize=TRUE)
raw_data_rma_bg_norm <- rma(rawData, background=TRUE, normalize=TRUE)

boxplot_raw <- boxplot(exprs(raw_data_rma), main="Boxplot raw data (Not normalized and no bg corr)")
boxplot_raw_bg <- boxplot(exprs(raw_data_rma_bg), main="Boxplot raw data (Not normalised but with bg corr")
boxplot_raw_norm <- boxplot(exprs(raw_data_rma_norm), main="Boxplot raw data (Normalized but no bg corr)")
boxplot_raw_bg_norm <- boxplot(exprs(raw_data_rma_bg_norm), main="Boxplot raw data (Normalized and with bg corr)")

plotMDS(raw_data_rma, main="MDS Plot raw data (Not normalized and without bg corr)")
plotMDS(raw_data_rma_bg, main="MDS Plot raw data (Not normalised but with bg corr)")
plotMDS(raw_data_rma_norm, main="MDS Plot raw data (Normalized but no bg corr)")
plotMDS(raw_data_rma_bg_norm, main="MDS Plot raw data (Normalized and with bg corr)")


### Annotation ###

head(fData(raw_data_rma_bg_norm)) # fData extracts feature data
rownames(head(fData(raw_data_rma_bg_norm)))

id <- featureNames(raw_data_rma_bg_norm)
symbol <- getSYMBOL(id, "hgu95av2.db")
fData(raw_data_rma_bg_norm) <- data.frame(symbol = symbol) #Ett sätt att lägga till en kolumn till en dataframe!
data_anno <- fData(raw_data_rma_bg_norm)

head(data_anno) #head shows first 6 elements

#Check which genes are AFFX genes and remove these as they are reference genes only
ref_names <- str_subset(row.names(data_anno), "^AFFX") #AFFX reference probes, list of names
ref_genes <- match(ref_names, row.names(data_anno)) #AFFX ref probes, list of indices (they are the last 67 elements!)
raw_data_rma_bg_norm_clean <- raw_data_rma_bg_norm[-ref_genes,] #CLEAN raw data!

### Data analysis! ###

design <- model.matrix(~0 + targets$estrogen)
colnames(design) <- c("absent", "present")
cont <- makeContrasts(contrasts = "present - absent", levels = design)

#Make linear fit for probes

fit <-lmFit(raw_data_rma_bg_norm_clean, design=design) #Make linear fit
fit <-contrasts.fit(fit, contrasts=cont) #Apply contrasts
fit <-eBayes(fit) #Calculate standard errors

fit_genes <- fit$genes$symbol #Char vector with all gene names

### Results!

#The topTable() function can be used to look at the top-ranked genes.

n_prob = nrow(fit) #number of probes (12558)
gene_table <- topTable(fit, number=nrow(fit)) #Make table with all probes in order

#Check number of significant genes, plus up/down
n <- topTable(fit, number=nrow(fit), p.value=0.05) # nrow(n) = 459 (sign genes)
dt_sign <- decideTests(fit) #A variable for which genes are up/down-regulated (1, 0 or -1)
table(dt_sign)

#dt_sign
#-1     0     1 
#176 12099   283 

#The 283 genes with the lowest adj p values are the down-regulated ones, 
#and the 176 genes with the highest values are up-regulated (after removing AFFX)

#Volcano plot:
volcanoplot(fit, highlight=459, names=fit_genes, main="Volcano plot for significant genes, p<0.05")














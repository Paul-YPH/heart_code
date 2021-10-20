###Processing of the data

library(Seurat)
library(dplyr)
list.files("F:/jky/jky10cluster/matrix/C2/")
C2.data<-Read10X(data.dir="F:/jky/jky10cluster/matrix/C2/")
C2 <- CreateSeuratObject(counts = C2.data,min.cells = 3, min.features = 200)
C2[["percent.mt"]] <- PercentageFeatureSet(C2, pattern = "^mt-")
VlnPlot(C2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(C2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(C2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
C2 <- subset(C2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
C2 <- NormalizeData(C2, normalization.method = "LogNormalize", scale.factor = 10000)
C2 <- FindVariableFeatures(C2, selection.method = "vst", nfeatures = 2000)


E2.data<-Read10X(data.dir="F:/jky/jky10cluster/matrix/E2/")
E2 <- CreateSeuratObject(counts = E2.data,min.cells = 3, min.features = 200)
E2[["percent.mt"]] <- PercentageFeatureSet(E2, pattern = "^mt-")
VlnPlot(E2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(E2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(E2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
E2 <- subset(E2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
#Normalizing
E2 <- NormalizeData(E2, normalization.method = "LogNormalize", scale.factor = 10000)
E2 <- FindVariableFeatures(E2, selection.method = "vst", nfeatures = 2000)


F2.data<-Read10X(data.dir="F:/jky/jky10cluster/matrix/F2/")
F2 <- CreateSeuratObject(counts = F2.data,min.cells = 3, min.features = 200)
F2[["percent.mt"]] <- PercentageFeatureSet(F2, pattern = "^mt-")
VlnPlot(F2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(F2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(F2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
F2 <- subset(F2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
F2 <- NormalizeData(F2, normalization.method = "LogNormalize", scale.factor = 10000)
#Identification of highly variable features (feature selection)
F2 <- FindVariableFeatures(F2, selection.method = "vst", nfeatures = 2000)


H1.data<-Read10X(data.dir="F:/jky/jky10cluster/matrix/H1/")
H1 <- CreateSeuratObject(counts = H1.data,min.cells = 3, min.features = 200)
H1[["percent.mt"]] <- PercentageFeatureSet(H1, pattern = "^mt-")
VlnPlot(H1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(H1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(H1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
H1 <- subset(H1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
H1 <- NormalizeData(H1, normalization.method = "LogNormalize", scale.factor = 10000)
H1 <- FindVariableFeatures(H1, selection.method = "vst", nfeatures = 2000)


C2$orig.ident<-"14M"
E2$orig.ident<-"7M"
F2$orig.ident<-"C"
H1$orig.ident<-"3M"

immune.anchors <- FindIntegrationAnchors(object.list = list(C2,E2,F2,H1), dims = 1:30,k.anchor =5,k.filter = 100 )
cm <- IntegrateData(anchorset = immune.anchors, dims = 1:30)

DefaultAssay(cm) <- "integrated"
# Run the standard workflow for visualization and clustering
cm <- ScaleData(cm, verbose = FALSE)
cm <- RunPCA(cm, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
cm <- RunTSNE(cm, reduction = "pca", dims = 1:20)
cm <- FindNeighbors(cm, reduction = "pca", dims = 1:20)
cm <- FindClusters(cm, resolution = 0.8)
###monocle
library(dplyr)
library(Seurat)
library(monocle3)
six.combined<-data
DefaultAssay(six.combined)<-"integrated"
data <- six.combined@assays$integrated@data
pd <-  six.combined@meta.data
#the metadata have many rubbish info,we delete it
#new_pd<-dplyr::select(pd,orig.ident,nCount_RNA,nFeature_RNA,percent.mt,seurat_clusters,CMT,all,T.ident,M.ident,type,hahaha,mm)
new_pd<-dplyr::select(pd,orig.ident,nCount_RNA,nFeature_RNA,percent.mt,seurat_clusters,type)
ad <-  as.data.frame(six.combined@active.ident)
try1<-cbind(new_pd,ad)
names(try1)[5] <- "cell_type"
head(try1)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
#create new cds obj
cds <- new_cell_data_set(data,cell_metadata  = try1,gene_metadata  = fData)
#save(cds,file = "cds.rdata")

cds <-preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds,preprocess_method = 'PCA',reduction_method = "UMAP",
                        umap.min_dist = 0.49,umap.n_neighbors = 40L,max_components =2)
cds<-cluster_cells(cds, reduction_method = "UMAP")
plot_cells(cds, label_groups_by_cluster=TRUE,  color_cells_by  = "six.combined@active.ident",cell_size = 0.5,group_label_size = 5)
cds <- learn_graph(cds,learn_graph_control = list(minimal_branch_len=5,euclidean_distance_ratio=0.7))
plot_cells(cds,  label_cell_groups=F, label_branch_points = F,label_leaves =F,label_roots =F,trajectory_graph_segment_size=0.7, color_cells_by = "mm",cell_size = 0.5,alpha = 1)+ scale_y_reverse()

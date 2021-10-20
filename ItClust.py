###Itclust
import ItClust as ic
import scanpy.api as sc
import os
from numpy.random import seed
from tensorflow import set_random_seed
import pandas as pd
import numpy as np
import warnings
os.environ["CUDA_VISIBLE_DEVICES"]="1"
warnings.filterwarnings("ignore")
#import sys
#!{sys.executable} -m pip install 'scanpy==1.4.4.post1'
#Set seeds
seed(20180806)
np.random.seed(10)
set_random_seed(20180806) # on GPU may be some other default

adata_train=sc.read("cmmm.h5ad") 
adata_test=sc.read("cmtmm.h5ad") 
adata_train.obs["celltype"]=adata_train.obs["type"] 
adata_test.obs["celltype"]=adata_train.obs["seurat_clusters"] 

clf=ic.transfer_learning_clf()  
clf.fit(adata_train, adata_test)
pred, prob, cell_type_pred=clf.predict()

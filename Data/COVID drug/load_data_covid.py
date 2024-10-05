"""
Created on Tue May  3 06:21:54 2022

@author: yrt05
"""

import numpy as np
import os.path as osp
import pyreadr
import pandas as pd
import torch
import os
import fsspec
# from sentence_transformers import SentenceTransformer
from torch_geometric.data import HeteroData, download_url, extract_zip
from torch_geometric.data import HeteroData
from torch_geometric.transforms import RandomLinkSplit, ToUndirected
from torch_geometric.nn import SAGEConv, to_hetero, GraphConv
from torch.nn import Linear
import torch.nn.functional as F
import torch.nn as nn
import torch_geometric.nn as pyg_nn
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from torch.nn import Parameter
from torch_geometric.nn import SGConv
from torch_geometric.nn import GAE, VGAE, GCNConv
from torch_geometric.nn import Sequential, HeteroConv
import random
random.seed(1)
torch.use_deterministic_algorithms(True)
np.random.seed(1)        
torch.manual_seed(1)
random_seed = 1
torch.manual_seed(random_seed)
torch.cuda.manual_seed(random_seed)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

#%% Functions

# ## "encoders=None": not encode the data
# # Just return "mapping"
def load_node_csv(path, index_col, encoders=None, **kwargs):   
    df = pd.read_csv(path, index_col=index_col, **kwargs)
    mapping = {index: i for i, index in enumerate(df.index.unique())}

    x = None
    if encoders is not None:
        xs = [encoder(df[col]) for col, encoder in encoders.items()]
        x = torch.cat(xs, dim=-1)

    return x, mapping


# ## "src_index_col": index columns of source nodes
# ## "dst_index_col": index columns of destination nodes
def load_edge_csv(path, src_index_col, src_mapping, dst_index_col, dst_mapping,
                  encoders=None, **kwargs):
    df = pd.read_csv(path, **kwargs)

    src = [src_mapping[index] for index in df[src_index_col]]
    dst = [dst_mapping[index] for index in df[dst_index_col]]
    edge_index = torch.tensor([src, dst])

    edge_attr = None
    if encoders is not None:
        edge_attrs = [encoder(df[col]) for col, encoder in encoders.items()]
        edge_attr = torch.cat(edge_attrs, dim=-1)

    return edge_index, edge_attr


#%% My data
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
# os.chdir("//COVID drug")
os.chdir("//COVID drug")

# Bind_data = pd.read_csv('BindingDB_All.tsv', sep='\t')



### test = pyreadr.read_r('DDI_dat_frame.Rdata')['DDI_dat_frame']


# DDI = (pyreadr.read_r('DDI_dat_frame.Rdata')['DDI_dat_frame']) # DDI   340 drugs (from 468 drugs)
# DTI = (pyreadr.read_r('DTI_dat_frame.Rdata')['DTI_dat_frame']) # DTI   214 d * 64 t
# PPI = (pyreadr.read_r('PPI_dat_frame.Rdata')['PPI_dat_frame']) # PPI   212 prot (from 495 genes)


DDI = (pyreadr.read_r('DDI_dat_frame.Rdata')['DDI_dat_frame']) # DDI   3410 drugs
DTI = (pyreadr.read_r('DTI_dat_frame.Rdata')['DTI_dat_frame']) # DTI   835 nodes, 825 edges




PPI = (pyreadr.read_r('PPI_dat_frame.Rdata')['PPI_dat_frame']) # PPI   212 prot (from 495 genes)



PPI_path = '//COVID drug//PPI_dat_frame.csv'
print(pd.read_csv(PPI_path).head())
PPI_x, PPI_mapping = load_node_csv(PPI_path, index_col='prot')

DDI_path = '//COVID drug//DDI_dat_frame.csv'
print(pd.read_csv(DDI_path).head())
DDI_x, DDI_mapping = load_node_csv(DDI_path, index_col='parent_key')

#%% Add indices
# =============================================================================
# for key in DDI_mapping:
#     DDI_mapping[key] += len(PPI_mapping) ### 0~962 : in total 963
# =============================================================================
#%%

PPI_feature = (  pyreadr.read_r("protein_seq_similarity.RData")  ['protein_seq_similarity']).values # PPI smilarity feature 659 x 659
DDI_feature = (  pyreadr.read_r("jaccard_sim_468_DDI.Rdata")  ['jaccard_sim_468_DDI']).values # DDI smilarity feature 3410 x 3410

drug_prot_RWR = (  pyreadr.read_r("drug_prot_RWR.Rdata")  ['drug_prot_RWR']).values # RWR = drug + prot



#%% PCA for Drug
x_drug = (  pyreadr.read_r("x_drug.Rdata")  ['x']).values 
print(x_drug.shape)

x_drug = np.concatenate((x_drug, DDI_feature), axis = 1)  # concat:  x_drug, DDI_feature
print(x_drug.shape)

# num_PCA = int(x_drug.shape[1]*2/3)
num_PCA_drug = 20

# num_take = 10
pca_1_drug = PCA(svd_solver = 'full',n_components = num_PCA_drug) # show first 10 PCs
pipe_1 = Pipeline([('scaler', StandardScaler()), ('pca', pca_1_drug)])

x_drug.shape

X = x_drug
# X = x_prot

X_transformed = pipe_1.fit_transform(X)

# X_transformed.shape

main_drug = pca_1_drug.explained_variance_ratio_



sum(main_drug[0:20])


x_drug = X_transformed



#%% PCA for Prot

PPI_feature.shape

x_prot = (  pyreadr.read_r("x_prot.RData")  ['prot_property']).values 
print(x_prot.shape)

x_prot = np.concatenate((x_prot, PPI_feature), axis = 1)  # concat:  x_prot, PPI_feature

print(x_prot.shape)

# num_PCA = int(x_prot.shape[1]*2/3)
# num_PCA = 350
num_PCA = 300

# num_take = 10
pca_1 = PCA(svd_solver = 'full',n_components = num_PCA) # show first 10 PCs
pipe_1 = Pipeline([('scaler', StandardScaler()), ('pca', pca_1)])
x_prot.shape
X = x_prot

X_transformed = pipe_1.fit_transform(X)

# X_transformed.shape

main = pca_1.explained_variance_ratio_

sum(main[0:300])

np.sum(main)
plt.style.use('default')


x_prot = X_transformed




#%% PCA plots
fig = plt.figure()
# fig.tight_layout() 
plt.subplots_adjust(hspace= 1)
plt.subplot(221)
plt.plot(main_drug,'o--', label = 'Explained Variance of each PC')
plt.xlabel('PCs', fontsize = 18)
plt.ylabel('Variance', fontsize = 18)
plt.grid()
plt.bar(range(0,num_PCA_drug), pca_1_drug.explained_variance_ratio_,
        alpha=0.5,
        align='center')
plt.step(range(0,num_PCA_drug), np.cumsum(pca_1_drug.explained_variance_ratio_),
         where='mid',
         color='red', label = 'Accumulative Explained Variance')
plt.title('(b) Drug Data', fontsize = 16, y=-1.5)
plt.legend(loc = 'best', bbox_to_anchor =(2, 2), fontsize = 16 )
# plt.show()





ax = plt.subplot(222)

plt.plot(main,'o--', label = 'Explained Variance of each PC')
# plt.axis(fontsize = 18)
plt.xlabel('PCs', fontsize = 18)
# plt.ylabel('Variance', fontsize = 18)
plt.grid()
plt.bar(range(0,num_PCA), pca_1.explained_variance_ratio_,
        alpha=0.5,
        align='center')
plt.step(range(0,num_PCA), np.cumsum(pca_1.explained_variance_ratio_),
         where='mid',
         color='red', label = 'Accumulative Explained Variance')
plt.title('(c) Protein Data', fontsize = 16, y=-1.5)
# plt.legend(loc = 'center right')

# ax.legend(bbox_to_anchor=(1.1, 1.05), fontsize = 16)
# fig.tight_layout(pad=5.0)
plt.show()




#%% Hete graph

class IdentityEncoder(object):
    # The 'IdentityEncoder' takes the raw column values and converts them to
    # PyTorch tensors.
    def __init__(self, dtype=None):
        self.dtype = dtype

    def __call__(self, df):
        return torch.from_numpy(df.values).view(-1, 1).to(self.dtype)


# All_path = 'C://Users//yrt05//Desktop//Covid_coference_presentation//edge_all.csv'
All_path = '//COVID drug//edge_all_label.csv'

def load_edge_csv_part(path, src_index_col, src_mapping, dst_index_col, dst_mapping,begin, end, 
                  encoders=None, **kwargs):
    df = pd.read_csv(path, **kwargs)[begin:end]

    src = [src_mapping[index] for index in df[src_index_col]]
    dst = [dst_mapping[index] for index in df[dst_index_col]]
    edge_index = torch.tensor([src, dst])

    edge_attr = None
    if encoders is not None:
        edge_attrs = [encoder(df[col]) for col, encoder in encoders.items()]
        edge_attr = torch.cat(edge_attrs, dim=-1)

    return edge_index, edge_attr


# all_edge_type = c("DDI", "PPI", "DTI")

edge_index_ddi, edge_label_ddi = load_edge_csv_part(
    All_path,
    src_index_col='V1',
    src_mapping=DDI_mapping,
    dst_index_col='V2',
    dst_mapping=DDI_mapping,
    begin = 0,
    end = 48460,
    # encoders={'label_type': IdentityEncoder(dtype=torch.long)}
    # encoders={'edge_type': IdentityEncoder(dtype=torch.long)},
)

edge_index_ppi, edge_label_ppi = load_edge_csv_part(
    All_path,
    src_index_col='V1',
    src_mapping=PPI_mapping,
    dst_index_col='V2',
    dst_mapping=PPI_mapping,
    begin = 48460,
    end = 51405, 
    # encoders={'label_type': IdentityEncoder(dtype=torch.long)},
    # encoders={'edge_type': IdentityEncoder(dtype=torch.long)},
)

edge_index_dti, edge_label_dti = load_edge_csv_part(
    All_path,
    src_index_col='V1',
    src_mapping=DDI_mapping,
    dst_index_col='V2',
    dst_mapping=PPI_mapping,
    begin = 51405,
    end = 51700,
    # encoders={'label_type': IdentityEncoder(dtype=torch.long)},
    # encoders={'edge_type': IdentityEncoder(dtype=torch.long)},
)


#%%



data = HeteroData()
data['Prot'].num_nodes = len(x_prot) 
data['Prot'].x = torch.from_numpy( x_prot).float()
data['Drug'].num_nodes = len(x_drug) 
data['Drug'].x = torch.from_numpy( x_drug).float()


data['Prot', 'interact', 'Prot'].edge_index = edge_index_ppi
# data['Prot', 'interact', 'Prot'].edge_label = edge_label_ppi

data['Drug', 'interact', 'Drug'].edge_index = edge_index_ddi

data['Drug', 'bind', 'Prot'].edge_index = edge_index_dti


# # edge_type_label


data_undirected = data

# data_undirected['Drug']
# data['Drug', 'bind', 'Prot']

import pickle
# Path to the pickle file
file_path_index = "data_06_20_index.pickle"
# file_path = "data_06_20.pickle"

# # Save the data to the pickle file
# with open(file_path, "wb") as file:
#     pickle.dump(data_undirected, file)

# Save the data to the pickle file
with open(file_path_index, "wb") as file:
    pickle.dump(data_undirected, file)
    
    
    
data = ToUndirected(merge=False)(data)

metadata = data.metadata()



train_data, val_data, test_data = RandomLinkSplit(is_undirected=True,edge_types=data_undirected.edge_types,)(data_undirected)

train_data, val_data, test_data = RandomLinkSplit(is_undirected=True,edge_types=data.edge_types,)(data)

train_data = ToUndirected(merge=True)(train_data)




# data_undirected












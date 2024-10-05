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

os.chdir("/Zika virus drug")


PPI_path = '/Zika virus drug//PPI_dat_frame.csv'
print(pd.read_csv(PPI_path).head())
PPI_x, PPI_mapping = load_node_csv(PPI_path, index_col='prot')


DDI_path = '/Zika virus drug//DDI_dat_frame.csv'
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


#%% PCA for Drug
x_drug = (  pyreadr.read_r("x_drug.Rdata")  ['x']).values 
print(x_drug.shape)

x_drug = np.concatenate((x_drug, DDI_feature), axis = 1)  # concat:  x_drug, DDI_feature
print(x_drug.shape)

# num_PCA_drug = int(x_drug.shape[1]*2/30)
num_PCA_drug = 50
# num_PCA_drug = 20

# num_take = 10
pca_1_drug = PCA(svd_solver = 'full',n_components = num_PCA_drug) # show first 10 PCs
pipe_1 = Pipeline([('scaler', StandardScaler()), ('pca', pca_1_drug)])

x_drug.shape

X = x_drug
# X = x_prot

X_transformed = pipe_1.fit_transform(X)

# X_transformed.shape

main_drug = pca_1_drug.explained_variance_ratio_




## Actual PCA for drug
drug_dim = 10


pca_drug = PCA(svd_solver = 'full',n_components = drug_dim) # show first 10 PCs
pipe_drug = Pipeline([('scaler', StandardScaler()), ('pca', pca_drug)])
x_drug_transformed = pipe_drug.fit_transform(x_drug)


x_drug = x_drug_transformed
sum(main_drug[0:drug_dim])


#%% PCA for Prot

PPI_feature.shape

x_prot = (  pyreadr.read_r("x_prot.RData")  ['prot_property']).values 
print(x_prot.shape)


x_prot  ### 343 prot amino acid
PPI_feature  ### 659 x 659 similarity

x_prot = np.concatenate((x_prot, PPI_feature), axis = 1)  # concat:  x_prot, PPI_feature

print(x_prot.shape)

num_PCA = min(x_prot.shape[0], int(x_prot.shape[1] *2/3) ) - 1

pca_1 = PCA(svd_solver = 'full',n_components = num_PCA) # show first 10 PCs
pipe_1 = Pipeline([('scaler', StandardScaler()), ('pca', pca_1)])
x_prot.shape
X = x_prot

X_transformed = pipe_1.fit_transform(X)

# X_transformed.shape

main_prot = pca_1.explained_variance_ratio_






## Actual PCA for prot
len(main_prot)
prot_dim = 550

pca_prot = PCA(svd_solver = 'full',n_components = prot_dim) # show first 10 PCs
pipe_prot = Pipeline([('scaler', StandardScaler()), ('pca', pca_prot)])
x_prot_transformed = pipe_prot.fit_transform(x_prot)


x_prot = x_prot_transformed
sum(main_prot[0:prot_dim])





#%% PCA plots


%matplotlib qt



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
plt.title('Drug Data', fontsize = 16)
plt.legend(loc = 'upper center', bbox_to_anchor =(1, -0.5), fontsize = 16 )

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)




ax = plt.subplot(222)

plt.plot(main_prot,'o--', label = 'Explained Variance of each PC')
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

plt.title('Protein Data', fontsize = 16)

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

plt.show()




#%% Edge 

All_path = '/Zika virus drug/edge_all_label.csv'

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
    end = 504534,
    # encoders={'label_type': IdentityEncoder(dtype=torch.long)}
    # encoders={'edge_type': IdentityEncoder(dtype=torch.long)},
)

edge_index_ppi, edge_label_ppi = load_edge_csv_part(
    All_path,
    src_index_col='V1',
    src_mapping=PPI_mapping,
    dst_index_col='V2',
    dst_mapping=PPI_mapping,
    begin = 504534  ,
    end = 2588+504534, 
    # encoders={'label_type': IdentityEncoder(dtype=torch.long)},
    # encoders={'edge_type': IdentityEncoder(dtype=torch.long)},
)

edge_index_dti, edge_label_dti = load_edge_csv_part(
    All_path,
    src_index_col='V1',
    src_mapping=DDI_mapping,
    dst_index_col='V2',
    dst_mapping=PPI_mapping,
    begin = 2588+504534 ,
    end = 507948  ,
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

data['Drug', 'interact', 'Drug'].edge_index = edge_index_ddi

data['Drug', 'bind', 'Prot'].edge_index = edge_index_dti


data_undirected = data

import pickle
# Path to the pickle file
file_path_index = "data_06_20_index.pickle"

# Save the data to the pickle file
with open(file_path_index, "wb") as file:
    pickle.dump(data_undirected, file)
    
    
    
data = ToUndirected(merge=False)(data)

metadata = data.metadata()


train_data, val_data, test_data = RandomLinkSplit(is_undirected=True,edge_types=data.edge_types,)(data)



























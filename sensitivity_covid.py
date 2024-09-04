# -*- coding: utf-8 -*-

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch
from torch_geometric.nn import HeteroConv, GraphConv, Linear, GCNConv
from torch.nn import Parameter
import pickle
import torch.nn as nn
from torch_geometric.utils import to_dense_adj
from torch_geometric.transforms import RandomLinkSplit, ToUndirected
from torch_geometric.utils import negative_sampling
from torch_geometric.utils import degree
from torch_geometric.transforms import NormalizeFeatures, Compose

### load


with open("Data/data_06_20_index_covid.pickle", 'rb') as f:
# with open("C:\\Users\\roy20001\\Desktop\\data_06_20_index.pickle", 'rb') as f:
    data_temp = pickle.load(f)

import os

os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'


device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
device = torch.device('cpu')


#%% To_undirect --> Split



data_temp = ToUndirected(merge=False)(data_temp) ### same edges
del data_temp['Prot', 'rev_interact', 'Prot'].edge_label  # Remove "reverse" label.
del data_temp['Drug', 'rev_interact', 'Drug'].edge_label  # Remove "reverse" label.
del data_temp['Prot', 'rev_bind', 'Drug'].edge_label  # Remove "reverse" label.

data_temp['Prot'].x

node_types = ['Prot', 'Drug']
transform_dict = {node_type: NormalizeFeatures() for node_type in node_types}
# transform = Compose(transform_dict)
# Apply normalization to each node type
for node_type in node_types:
    node_features = data_temp[node_type].x
    normalized_features = (node_features - node_features.mean(dim=0)) / node_features.std(dim=0)
    data_temp[node_type].x = normalized_features


data_temp['Prot'].x


train_data, val_data, test_data = RandomLinkSplit(is_undirected=False,
                                                  edge_types=data_temp.edge_types,
                                                  )(data_temp)



metadata = train_data.metadata()

x_prot = train_data['Prot']
x_drug = train_data['Drug']

group_name_concat = []
for item in metadata[1]:
    string = '__'.join(item)
    group_name_concat.append(string)




x_dict = train_data.x_dict

edge_index_dict = train_data.edge_index_dict
edge_label_index = train_data['Drug', 'Prot'].edge_label_index
# train_data['Prot', 'Drug'].edge_label_index


# global adj_mat_sub
data = train_data
edge_index_dict = data.edge_index_dict
import os


os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

# device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

seed = 333
torch.manual_seed(seed)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(seed)
import random
random.seed(seed)
# import numpy as np
import numpy as np
np.random.seed(seed)


# %% To_undirect --> Split

data_temp = ToUndirected(merge=False)(data_temp)  # same edges
# Remove "reverse" label.
del data_temp['Prot', 'rev_interact', 'Prot'].edge_label
# Remove "reverse" label.
del data_temp['Drug', 'rev_interact', 'Drug'].edge_label
del data_temp['Prot', 'rev_bind', 'Drug'].edge_label  # Remove "reverse" label.

data_temp['Prot'].x

node_types = ['Prot', 'Drug']
transform_dict = {node_type: NormalizeFeatures() for node_type in node_types}
# transform = Compose(transform_dict)



#%% Load model
def find_Adj(edge_index_dict):

    # edge_index_dict = data.edge_index_dict

    # PPI
    edge_PPI = edge_index_dict['Prot', 'interact', 'Prot']
    A_PPI = to_dense_adj(edge_PPI, max_num_nodes=x_prot['num_nodes'])[
        0].clamp(max=1)
    edge_PPI_t = edge_index_dict['Prot', 'rev_interact', 'Prot']
    A_PPI_t = to_dense_adj(edge_PPI_t, max_num_nodes=x_prot['num_nodes'])[
        0].clamp(max=1)
    A_PPI = (A_PPI + A_PPI_t).clamp(max=1)

    # DDI
    edge_DDI = edge_index_dict['Drug', 'interact', 'Drug']
    A_DDI = to_dense_adj(edge_DDI, max_num_nodes=x_drug['num_nodes'])[
        0].clamp(max=1)
    edge_DDI_t = edge_index_dict['Drug', 'rev_interact', 'Drug']
    A_DDI_t = to_dense_adj(edge_DDI_t, max_num_nodes=x_drug['num_nodes'])[
        0].clamp(max=1)
    A_DDI = (A_DDI + A_DDI_t).clamp(max=1)

    # Self loop
    A_PPI = (A_PPI + torch.eye(x_prot['num_nodes'])).clamp(max=1)
    A_DDI = (A_DDI + torch.eye(x_drug['num_nodes'])).clamp(max=1)

    # DTI
    edge_pair_DTI = edge_index_dict['Drug',
                                    'bind', 'Prot']  # DTI : drug [0], prot [1]
    edge_pair_DTI[0] = edge_pair_DTI[0] + x_prot['num_nodes']
    adj_mat_sub_DTI = to_dense_adj(
        edge_pair_DTI, max_num_nodes=x_prot['num_nodes'] + x_drug['num_nodes'])[0].clamp(max=1)
    edge_pair_DTI[0] = edge_pair_DTI[0] - x_prot['num_nodes']

    # IDI
    # DTI_t :  prot [0], drug [1]
    edge_pair_DTI_t = edge_index_dict['Prot', 'rev_bind', 'Drug']
    edge_pair_DTI_t[1] = edge_pair_DTI_t[1] + x_prot['num_nodes']
    adj_mat_sub_DTI_t = to_dense_adj(
        edge_pair_DTI_t, max_num_nodes=x_prot['num_nodes'] + x_drug['num_nodes'])[0].clamp(max=1)
    edge_pair_DTI_t[1] = edge_pair_DTI_t[1] - x_prot['num_nodes']

    Adj_all = torch.zeros(adj_mat_sub_DTI_t.shape)

    Adj_all[0:(A_PPI.shape[0]), 0:(A_PPI.shape[1])] = A_PPI
    Adj_all[(A_PPI.shape[0]):, (A_PPI.shape[1]):] = A_DDI

    Adj_all[(A_PPI.shape[1]):, 0:(A_PPI.shape[0])
            ] = adj_mat_sub_DTI[(A_PPI.shape[1]):, 0:(A_PPI.shape[0])]
    Adj_all[0:(A_PPI.shape[0]), (A_PPI.shape[1])
               :] = adj_mat_sub_DTI_t[0:(A_PPI.shape[0]), (A_PPI.shape[1]):]

    return Adj_all

class Identity(nn.Module):
    def forward(self, x):
        return x
class GraphConvolution(nn.Module):
    def __init__(self, hidden_channels, out_features, num_layer):

        super(GraphConvolution, self).__init__()

        random_seed = 321

        torch.manual_seed(random_seed)
        torch.cuda.manual_seed(random_seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False

        in_features = hidden_channels
        out_features = hidden_channels

        num_layers_before_merging = 1
        num_layers_after_merging = num_layer - num_layers_before_merging

        self.leaky_relu = nn.LeakyReLU()
        # %% Before merging
        self.encoder_before_merging = torch.nn.ModuleList()
        for _ in range(num_layers_before_merging):
            conv = HeteroConv({
                edge_type: GraphConv(-1, hidden_channels,
                                     bias = False)  # GraphConv
                # edge_type: GCNConv(-1, hidden_channels, cached = True, add_self_loops = False)  ### GCNConv
                for edge_type in metadata[1]
            })
            self.encoder_before_merging.append(conv)  # Add function here

        # %% After merging

        # self.hidden_layers_after_merge = nn.ModuleList([GCN(hidden_channels, hidden_channels) for _ in range(num_layers_after_merging)])
        self.num_layers_after_merging = num_layers_after_merging
        self.weight_matrices = nn.ParameterList()
        for i in range(num_layers_before_merging + num_layers_after_merging):
            if i == 0:
                self.weight_matrices.append(nn.Parameter(
                    torch.Tensor(in_features, hidden_channels)))
            elif i != num_layers_before_merging + num_layers_after_merging-1:
                self.weight_matrices.append(nn.Parameter(torch.Tensor(
                    hidden_channels, hidden_channels)))  # middle layer
            else:
                self.weight_matrices.append(nn.Parameter(
                    torch.Tensor(hidden_channels, out_features)))  # last layer


        self.layer_norm_prot = Identity()
        self.layer_norm_drug = Identity()

        self.layer_norm_identity = Identity()
        # self.layer_norm = nn.BatchNorm2d(hidden_channels)
        # self.layer_norm_last = Identity()

        self.reset_parameters()
        self.intermediate_representations = []

    def reset_parameters(self):
        for weight in self.weight_matrices:


            nn.init.xavier_uniform_(weight)

# %% Forward path

    def forward(self, z_dict, edge_index_dict, edge_label_index):
        # def forward(self, z_dict, edge_index_dict):
        # %% First layer (before merging)

        # z_dict = x_dict
        z_dict_normalized = z_dict
        # z_dict_normalized = {key: self.layer_norm(x) for key, x in z_dict_normalized.items()}
        z_dict_normalized['Drug'] = self.layer_norm_drug(z_dict_normalized['Drug'])
        z_dict_normalized['Prot'] = self.layer_norm_prot(z_dict_normalized['Prot'])


        self.intermediate_representations = []  # latent features

        # %% Merging before (for same dimension)
        for conv in self.encoder_before_merging:  # encoder_before_merging = 1 for same dimension
            z_dict_normalized = conv(z_dict_normalized, edge_index_dict)


        # After merging
        # Feature: stack (drug + prot)
        z = torch.cat([z_dict_normalized['Prot'],
                      z_dict_normalized['Drug']], dim=0)

        # record first (merging) layer
        self.intermediate_representations.append(z)  # H_0

        # Large Adj matrix (with projection)

        # Graph feature
        graph_x_prot = z[0: (z_dict_normalized['Prot'].shape[0]), :]
        graph_x_drug = z[(z_dict_normalized['Prot'].shape[0])                         :, :]  # graph prot features


        # No graph features
        graph_x_prot = torch.zeros(graph_x_prot.shape)
        graph_x_drug = torch.zeros(graph_x_drug.shape)

        graph_property = torch.cat((graph_x_prot, graph_x_drug), dim=0)

        # %% A * X * W
        Adj_all = find_Adj(edge_index_dict)
        # Adj_all = Adj_all + torch.eye(Adj_all.shape[0], Adj_all.shape[1])

        # A normalize (mean aggregate)
        degrees_A_inv = torch.diag(torch.sum(Adj_all, dim=1)**(-1))
        Adj_all = degrees_A_inv @ Adj_all

        ### A * X
        z = Adj_all @ z  # torch.Size([963, 10]) project to same dimension

        all_x = torch.cat([z, graph_property], dim=0)

        # %% First layer (to put all_x together)
        ### AX * W

        all_x = self.layer_norm_identity(all_x @ self.weight_matrices[0])

        self.intermediate_representations.append(
            self.weight_matrices[0])  # W_0

        z_update = all_x[0: (z_dict_normalized['Prot'].shape[0] +
                             z_dict_normalized['Drug'].shape[0]), :]  # node features
        graph_x_update = all_x[(z_dict_normalized['Prot'].shape[0] +
                                z_dict_normalized['Drug'].shape[0]):, :]  # Graph feature

        self.intermediate_representations.append(z_update)  # H_1

        # torch.max(z_update) ### can have value greater than 1
        # z_update.shape

        # %% Following layers

        for i in range(0, self.num_layers_after_merging):

            z_node_all = Adj_all @ z_update
            # node + graph


            all_x = torch.cat([z_node_all, graph_x_update], dim=0)
            # all_x = (all_x @ weight_matrices[i+1]).relu()


            # print('layer number')

            if i == self.num_layers_after_merging-1:
            # if i >= self.num_layers_after_merging-2:

                all_x = (all_x @ self.weight_matrices[i+1])  ### +1 due to 1 layer projection


            else:
                all_x = all_x @ self.weight_matrices[i+1]
                all_x = (all_x).relu()


            # node features
            z_update = all_x[0: (z_dict['Prot'].shape[0] +
                                 z_dict['Drug'].shape[0]), :]
            # Graph feature
            graph_x_update = all_x[( z_dict['Prot'].shape[0] + z_dict['Drug'].shape[0] ):, :]

            # print('graph_x_update', graph_x_update)
            # record updated latent features
            # W_1,W_2,...
            self.intermediate_representations.append(self.weight_matrices[i+1])
            self.intermediate_representations.append(z_update)  # H_2,H_2,...

        return z_update


class EdgeDecoder(torch.nn.Module):
    def __init__(self, hidden_channels):
        super().__init__()

    def forward(self, z, edge_label_index):
        # z = z_encoded
        z_update = z  # Prot + Drug   z_update.shape  z = z

        output = torch.sigmoid(torch.matmul(z_update, z_update.t()))

        A_adj_predicted = output

        return A_adj_predicted


# %% Model integration

class Model(torch.nn.Module):
    def __init__(self, hidden_channels, num_layers):
        super().__init__()
        self.encoder = GraphConvolution(
            hidden_channels, hidden_channels, num_layers)
        # self.encoder = to_hetero(self.encoder, data.metadata(), aggr='mean')
        self.decoder = EdgeDecoder(hidden_channels)

    def forward(self, x_dict, edge_index_dict, edge_label_index):  # "edge_label_index" not used
        z_dict = self.encoder(x_dict, edge_index_dict, edge_label_index)
        return self.decoder(z_dict, edge_label_index),   self.encoder.intermediate_representations
        # return self.decoder(z_dict, edge_label_index)

    # decompose "forward" into two parts

    def encode(self, x_dict, edge_index_dict, edge_label_index):
        z_dict = self.encoder(x_dict, edge_index_dict, edge_label_index)
        return z_dict

    def decode(self, z, edge_label_index_combine, edge_label_combine):  # "edge_label_combine" not used
        # same as inner product
        specific_edge = ((z[edge_label_index_combine[0]] *
                         z[edge_label_index_combine[1]]).sum(dim=-1))
        specific_edge = torch.sigmoid(specific_edge)
        return specific_edge

### model from 'model.pth'

# Load the saved model and data
file_path = 'model_covid.pth'
checkpoint = torch.load(file_path)

saved_state_dict = torch.load(file_path)


# num_node = 3410 + 659
num_layers = 5
num_layer = num_layers

num_layers_before_merging = 1
# num_layers_after_merging = num_layers - num_layers_before_merging
# in_features = 963 + 2
hidden_channels = 16

# Modify the keys in the saved state dictionary to include angle brackets
modified_state_dict = {}
for key, value in saved_state_dict.items():
    # Modify the keys to include angle brackets
    modified_key = key.replace('__', '<').replace('_', '>')

    # Store the modified key-value pair
    modified_state_dict[modified_key] = value

loaded_model = Model(hidden_channels=hidden_channels, num_layers=num_layers).to(device)

for key in modified_state_dict:
    new_key = key.replace('<', '__').replace('>', '_')
    loaded_model.state_dict()[new_key] = modified_state_dict[key]

data_in = data_temp

pred, _ = loaded_model(data_in.x_dict, data_in.edge_index_dict, 0)

#%% My net

H_list = []

indices_to_add = []

#print(len(_))

for i in range(num_layer + 1):
  indices_to_add.append(2*i)
#indices_to_add = [0, 2, 4, 6, 8]
for index in indices_to_add:
    H_list.append(_[index])
H_list

### For W_list
W_list = []
indices_to_add = []

#indices_to_add = [1,3,5,7]
for i in range(num_layer):
  indices_to_add.append(2*i + 1)

for index in indices_to_add:
    W_list.append(_[index])
W_list
# _[9]

ground_truth = find_Adj(data_in.edge_index_dict)

def sparse_identity(n):
    """Create a sparse identity matrix of size n x n."""
    indices = torch.arange(n)
    indices = torch.stack((indices, indices), dim=0)
    values = torch.ones(n)
    return torch.sparse_coo_tensor(indices, values, (n, n))

def partial_log_term_partial_A_bar(H_li, H_lj, A_bar, A_hat_ij, mask, i, j):  ### equation (3)

    I_n = sparse_identity(len(A_bar))  # Identity matrix of size n


    #%% unsolved
    term1 = H_li_to_A_bar(A_bar, H_list, W_list, I_n, mask, i) @ kronecker(I_n, H_lj.unsqueeze(0).T.to_sparse())


    term2 = kronecker(I_n, H_li.unsqueeze(0).to_sparse()) @ H_lj_T_to_A_bar(A_bar, H_list, W_list, I_n, mask, j)

    return (1 - A_hat_ij) * (term1 + term2)

def safe_log(x, epsilon=1e-10):
    return torch.log(x + epsilon)


    #%% Eqn (4) : partial_Hi_to_A_bar
    ### transpose : partial_Hj_T_to_A_bar
def H_li_to_A_bar(A_bar, H_list, W_list, I_n, mask, i):

# =============================================================================
#         partial_A_i_to_A_bar = 0
# =============================================================================

    selector_ei = torch.zeros(A_bar.shape[0])
    selector_ei[i] = 1
    selector_ei = selector_ei.unsqueeze(0).to_sparse().T

    n = A_bar.shape[0]

    U_bar_nn = PermutationRelatedMatrix(n, n)

    partial_A_i_to_A_bar = kronecker( I_n, selector_ei.t() ) @ U_bar_nn

    l = num_layer

    part_1 = partial_A_i_to_A_bar @ ( kronecker( I_n, (H_list[ l -1 ] @ W_list[l-1]).to_sparse()  ) ) ## H_1 W_2   H_0 == X_0

    #%% equation (4) part 2
    # H_l_minus_1_to_A_bar(A_bar, H_list, W_list, mask, I_n, l) ### 81 x 18
    part_2 = kronecker( I_n, A_bar[i,:].unsqueeze(0).to_sparse() ) @ H_l_minus_1_to_A_bar(A_bar, H_list, W_list, mask, I_n, l) @ kronecker( I_n, W_list[l-1].to_sparse() )


    braket = part_1 + part_2

    mask_current = mask[l-1  -1 ]

# =============================================================================
#         braket.shape[0]
#
#         mask_current.shape[0]
# =============================================================================

    # J_nn = torch.ones( int(braket.shape[0] / mask_current.shape[0]), int(braket.shape[1] / mask_current.shape[1])).to_sparse()

    J_nn = torch.ones(1,n).to_sparse()

    result =  sparse_sparse_hadamard_product(  kronecker(J_nn, mask_current.to_sparse()) ,   braket)

# =============================================================================
#         # result =  mask_current.to_sparse() @ braket
# =============================================================================

    return result

#%% dual H_lj_T_to_A_bar
def H_lj_T_to_A_bar(A_bar, H_list, W_list, I_n, mask, j):

    selector_ej = torch.zeros(A_bar.shape[0])
    selector_ej[j] = 1
    selector_ej = selector_ej.unsqueeze(0).to_sparse()
    selector_ej = selector_ej.T

    n = A_bar.shape[0]
    U_bar_nn = PermutationRelatedMatrix(n, n)

    l = num_layer

    # H_list[l][i,:].shape


    #%% part 1
    ### partial_H_l_minus_1_T_to_A_bar    18 x 81

    part_1 = kronecker( I_n, W_list[l-1].T.to_sparse() ) @ partial_H_l_minus_1_T_to_A_bar(A_bar, H_list, W_list, mask, I_n, l) @ kronecker( I_n, A_bar[j,:].unsqueeze(0).T.to_sparse() )

    #%% part 2

    part_2 = kronecker( I_n, ( W_list[l-1].T @ H_list[ l -1 ].T ).to_sparse() ) @ U_bar_nn @ kronecker( I_n , selector_ej )  ## H_1 W_2   H_0 == X_0

    braket = part_1 + part_2

    mask_current = mask[l-1  -1 ]

    J_nn = torch.ones(n,1).to_sparse()

    # J_nn = torch.ones( int(braket.shape[0] / mask_current.shape[0]), int(braket.shape[1] / mask_current.shape[1])).to_sparse()

    result =  sparse_sparse_hadamard_product(  kronecker(J_nn, mask_current.T.to_sparse() )  ,   braket)

# =============================================================================
#         # result =  mask_current.to_sparse() @ braket
# =============================================================================

    return result


#%% mask activation (layer 0 ~ l)


def mask_activation(A_bar, H_list, W_list, num_layer):

    mask_list = []

    for i in range(num_layer):
        #print("i: ", i)
        mask_list.append( torch.ge(A_bar @ H_list[i] @ W_list[i], 0).int() ) # Use .int() to convert boolean values to integer (1 or 0)
    return mask_list




#%% Iterative process involved
l = num_layer

# =============================================================================
# test = H_l_minus_1_to_A_bar(A_bar, H_list, W_list, mask, I_n, l)
# =============================================================================

def H_l_minus_1_to_A_bar(A_bar, H_list, W_list, mask, I_n, l):

# =============================================================================
#         mask_current = mask[l-1 -1] # last layer (l-1) never used
# =============================================================================

    num_layer = l
# =============================================================================
#         layer_i = num_layer - 1
# =============================================================================

# =============================================================================
#         mask_current = mask[l -1 -1]
# =============================================================================

    n = A_bar.shape[0]
    U_bar_nn = PermutationRelatedMatrix(n, n)



    #%% i = 1 , 2

    i = l-1  ### start from last layer

    part_1 = U_bar_nn @ kronecker( I_n, (H_list[ l-1 - i] @ W_list[ l-1 - i]).to_sparse())

    # partial_H_to_A = H_l_i_to_A_bar(A_bar, H_list, W_list, mask, I_n, layer_i, l) ### next iteration

    if num_layer - i == 1:
        partial_H_to_A = U_bar_nn @ kronecker( I_n, (H_list[1] @ W_list[0]).to_sparse())

    part_2 = kronecker( I_n, A_bar.to_sparse()) @ partial_H_to_A @ kronecker( I_n, W_list[ l-1 - i].to_sparse() )

    J_nn = torch.ones(n,n).to_sparse()

    # J_nn = torch.ones( n , W_list[ l-1 - i].shape[1]).to_sparse()

    ### mask
    mask_current = mask[ l-1 - i]
    # print('Mask: ', l-1 - i)
    # J_nn = torch.ones(n,n).to_sparse()

    result = sparse_sparse_hadamard_product( kronecker(J_nn, mask_current.to_sparse()) , (  (part_1 + part_2)  ) )

    for i in range(l-2, 0, -1):

        part_1 = U_bar_nn @ kronecker( I_n, (H_list[ l-1 - i] @ W_list[ l-1 - i]).to_sparse())

        # layer_i = i

        # partial_H_to_A = H_l_i_to_A_bar(A_bar, H_list, W_list, mask, I_n, layer_i, l) ### next iteration

        if num_layer - i == 1:
            partial_H_to_A = U_bar_nn @ kronecker( I_n, (H_list[1] @ W_list[0]).to_sparse())
        else:
            partial_H_to_A = result ### complete the loop

        part_2 = kronecker( I_n, A_bar.to_sparse()) @ partial_H_to_A @ kronecker( I_n, W_list[ l-1 - i].to_sparse() )

        # J_nn = torch.ones( n , W_list[ l-1 - i].shape[1]).to_sparse()

        ### mask
        J_nn = torch.ones(n,n).to_sparse()
        mask_current = mask[ l-1 - i]
        # print('Mask: ', l-1 - i)
        # J_nn = torch.ones(n,n).to_sparse()

        result = sparse_sparse_hadamard_product( kronecker(J_nn, mask_current.to_sparse()) ,  (part_1 + part_2)   )

    return result


#%% partial_H_l_minus_1_T_to_A_bar

# =============================================================================
# test = partial_H_l_minus_1_T_to_A_bar(A_bar, H_list, W_list, mask, I_n, l)
# =============================================================================

def partial_H_l_minus_1_T_to_A_bar(A_bar, H_list, W_list, mask, I_n, l):

    # num_layer = l

    n = A_bar.shape[0]
    U_nn = PermutationMatrix(n, n)
    #%%

    i = l-1  ### start from last layer

    if l - i == 1:

        partial_H_T_to_A = kronecker( I_n, (H_list[1] @ W_list[0]).T.to_sparse()) @ U_nn

    part_1 = kronecker( I_n, W_list[ l-1 - i].T.to_sparse() ) @ partial_H_T_to_A @ kronecker( I_n, A_bar.T.to_sparse())


    part_2 = kronecker( I_n, (H_list[ l-1 - i] @ W_list[ l-1 - i]).T.to_sparse()) @ U_nn


    ### mask
    mask_current = mask[ l-1 - i]

# =============================================================================
#         # J_nn = torch.ones( int(part_2.shape[0] / mask_current.shape[0]), int(part_2.shape[1] / mask_current.shape[1] )).to_sparse()
#         J_nn = torch.ones( n * W_list[ l-1 - i].T.shape[1], n * n).to_sparse()
#
# =============================================================================
    J_nn = torch.ones(n,n).to_sparse()


    # J_shape =  (  (part_1 + part_2) @ kronecker( I_n, W_list[l-1 - i].T.to_sparse() )  ).shape
    # J_nn = torch.ones(int ( J_shape[0] / mask_current.shape[0]), int ( J_shape[1] / mask_current.shape[1])).to_sparse()

    ### 81 x 18
    result = sparse_sparse_hadamard_product(  kronecker(J_nn, mask_current.to_sparse()).T  ,   (part_1 + part_2)   )

    for i in range(l-2, 0, -1):


        # i = 2

        if l - i == 1:

            partial_H_T_to_A = kronecker( I_n, (H_list[1] @ W_list[0]).T.to_sparse()) @ U_nn

        part_1 = kronecker( I_n, W_list[ l-1 - i].T.to_sparse() ) @ partial_H_T_to_A @ kronecker( I_n, A_bar.T.to_sparse())


        part_2 = kronecker( I_n, (H_list[ l-1 - i] @ W_list[ l-1 - i]).T.to_sparse()) @ U_nn


        ### mask
        mask_current = mask[ l-1 - i]


# =============================================================================
#             J_nn = torch.ones( n * W_list[ l-1 - i].T.shape[1], n * n).to_sparse()
# =============================================================================

        J_nn = torch.ones(n,n).to_sparse()

        # J_shape =  (  (part_1 + part_2) @ kronecker( I_n, W_list[l-1 - i].T.to_sparse() )  ).shape
        # J_nn = torch.ones(int ( J_shape[0] / mask_current.shape[0]), int ( J_shape[1] / mask_current.shape[1])).to_sparse()

        ### 81 x 18
        result = sparse_sparse_hadamard_product(  kronecker(J_nn, mask_current.to_sparse()).T  ,   (part_1 + part_2)   )

    return result





#%% kronecker : Sparse tensor
def kronecker(sparse_tensor1, sparse_tensor2):
    # Extract indices and values from the first sparse tensor
    indices1 = sparse_tensor1.coalesce().indices().t()
    values1 = sparse_tensor1.coalesce().values()
    indices2 = sparse_tensor2.coalesce().indices().t()
    values2 = sparse_tensor2.coalesce().values()

    result_indices = []
    result_values = []
    dim1 = sparse_tensor1.size()
    dim2 = sparse_tensor2.size()

    for idx1, v1 in zip(indices1, values1):
        for idx2, v2 in zip(indices2, values2):
            new_index = (idx1[0] * dim2[0] + idx2[0], idx1[1] * dim2[1] + idx2[1])
            new_value = v1 * v2
            result_indices.append(new_index)
            result_values.append(new_value)

    result_indices = torch.tensor(result_indices).t()
    result_values = torch.tensor(result_values)
    result_size = (dim1[0] * dim2[0], dim1[1] * dim2[1])
    return torch.sparse_coo_tensor(result_indices, result_values, result_size)

def sparse_sparse_hadamard_product(sparse1, sparse2):
    # Ensure both tensors are in COO format
    assert sparse1.is_sparse and sparse2.is_sparse

    # Find matching indices in both tensors
    indices1 = sparse1.coalesce().indices().t()
    values1 = sparse1.coalesce().values()
    indices2 = sparse2.coalesce().indices().t()
    values2 = sparse2.coalesce().values()

    # Dictionary to hold the results
    result_dict = {}

    # Convert indices to tuples to use them as dictionary keys
    for idx1, val1 in zip(indices1, values1):
        result_dict[tuple(idx1.numpy())] = val1

    # Perform element-wise multiplication for matching indices
    for idx2, val2 in zip(indices2, values2):
        idx_tuple = tuple(idx2.numpy())
        if idx_tuple in result_dict:
            result_dict[idx_tuple] *= val2
        else:
            result_dict[idx_tuple] = 0  # or omit this to maintain sparsity

    # Rebuild the sparse tensor from the result dictionary
    if result_dict:
        new_indices = torch.tensor(list(result_dict.keys())).t()
        new_values = torch.tensor(list(result_dict.values()))
        result_size = sparse1.size()  # Assuming both inputs are the same size
        return torch.sparse_coo_tensor(new_indices, new_values, result_size)
    else:
        return torch.sparse_coo_tensor(size=sparse1.size())




# %% Functions (sparse tensor)

def unitvector(dimension, position):
    """ Create a sparse tensor with a single non-zero element at the specified position. """
    indices = torch.tensor([[position], [0]], dtype=torch.int64)  # Position adjusted to 0-indexing
    values = torch.tensor([1.], dtype=torch.float32)
    shape = torch.Size([dimension, 1])
    return torch.sparse_coo_tensor(indices, values, shape)

def ElementaryMatrix(row, col, position1, position2):
    """ Create the elementary matrix as a sparse tensor. """
    u = unitvector(row, position1 - 1)  # Adjust for 0-indexing
    v = unitvector(col, position2 - 1).t()  # Transpose to make it a row vector
    # Multiplying a sparse and dense matrix, which is v converted to dense
    E = torch.sparse.mm(u, v.to_dense())
    return E

def PermutationMatrix(row, col):
    # Start with an empty sparse tensor
    indices = torch.empty([2, 0], dtype=torch.int64)
    values = torch.empty([0], dtype=torch.float32)
    shape = torch.Size([row * col, row * col])
    U = torch.sparse_coo_tensor(indices, values, shape)

    for i in range(row):
        for j in range(col):
            E_ij = ElementaryMatrix(row, col, i + 1, j + 1).to_dense()  # Conversion is necessary for kron
            E_ji = ElementaryMatrix(col, row, j + 1, i + 1).to_dense()
            # Kronecker product needs to be handled in dense due to current PyTorch limitations
            U += torch.kron(E_ij, E_ji).to_sparse()  # Convert result back to sparse

    return U

def PermutationRelatedMatrix(row, col):
    # Start with an empty sparse tensor
    indices = torch.empty([2, 0], dtype=torch.int64)
    values = torch.empty([0], dtype=torch.float32)
    shape = torch.Size([row * row, col * col])
    U = torch.sparse_coo_tensor(indices, values, shape)

    for i in range(row):
        for j in range(col):
            E_ij = ElementaryMatrix(row, col, i + 1, j + 1).to_dense()  # Conversion to dense
            # Kronecker product needs to be handled in dense
            U += torch.kron(E_ij, E_ij).to_sparse()  # Convert result back to sparse

    return U

# =============================================================================
# # Example usage
# row, col = 3, 3
# U_permutation = PermutationMatrix(row, col)
# U_permutation_related = PermutationRelatedMatrix(row, col)
# =============================================================================

#%%
#%% Derivative of the L to A_bar
### data from "#%% get H_list and W_list"
A_pred = pred

### A normalize (mean aggregate)
Adj_all = find_Adj(data_temp.edge_index_dict)
degree_node = torch.sum(Adj_all, dim=1)



#degrees_A_inv = torch.diag(torch.sum(Adj_all, dim=1)**(-1))
#A_bar = degrees_A_inv @ Adj_all
# A_bar = A

data_temp.metadata()

#%% Removal of nodes (indices)

#A_test = torch.rand(3,3)
A_test = find_Adj(data_temp.edge_index_dict)
Adj_all = A_test


A_no_self = Adj_all - torch.eye(Adj_all.size()[0])

degree_node = torch.sum(A_no_self, dim=1)

#node_degree_select = torch.sort(degree_node, descending = True)   ### highest degrees: default ascending sort
node_degree_select = torch.sort(degree_node, descending = False)   ### lowest degrees: default ascending sort

### node indices
node_indices = node_degree_select[1]
#print("Num nodes: ", len(node_indices))


# # Step 1: Mask out non-positive elements
# masked_degree = degree_node[degree_node > 0]

# Step 1: Create a mask for nodes with degree greater than 0
mask_positive_degree = node_degree_select[0] > 0

# Step 2: Create a mask for node indices smaller than 495 (total 963)
mask_smaller_indices = node_degree_select[1] < data_temp['Prot'].x.shape[0]

# Step 3: Combine the masks using logical AND
combined_mask = mask_positive_degree & mask_smaller_indices

# Step 4: Apply the combined mask to select node degrees and indices
selected_degrees = node_degree_select[0][combined_mask]
selected_indices = node_degree_select[1][combined_mask]

print('Number of Prot for Selection: ', len(selected_indices))

# Step 5: Find the n minimal node degrees and their corresponding indices
n = 100  # Number of minimal nodes to find
min_values, min_indices = selected_degrees.topk(k=n, largest=False)

# Select the corresponding indices from the selected indices

print('How many iters (10): ', len(selected_indices)/10)

iters = 20

for ind in range(0, iters):
# i = 0
    minimal_node_indices = selected_indices[min_indices[10*ind: 10*(ind+1)]]
    # minimal_node_indices = selected_indices[min_indices[10:20]]
# =============================================================================
#     ind = 1
# =============================================================================
    
    node_degree_select_indices = minimal_node_indices
    
    ### updated A matrix after node removal
    
    H_list = []
    
    indices_to_add = []
    
    #print(len(_))
    
    for i in range(num_layer + 1):
      indices_to_add.append(2*i)
    #indices_to_add = [0, 2, 4, 6, 8]
    for index in indices_to_add:
        H_list.append(_[index])
        #print(_[index])
    H_list
    
    ### For W_list
    W_list = []
    indices_to_add = []
    
    #indices_to_add = [1,3,5,7]
    for i in range(num_layer):
      indices_to_add.append(2*i + 1)
    
    for index in indices_to_add:
        W_list.append(_[index])
    W_list
    
    #%% select from H and W
    ## H_list = [a[0:few_index_test] for a in H_list]
    H_list = [a[node_degree_select_indices] for a in H_list]  ### select rows for each H in H_list
    
    #print("H_list", H_list)
    
    #W_list
    
    degrees_A_inv = torch.diag(torch.sum(Adj_all, dim=1)**(-1))
    A_bar = degrees_A_inv @ Adj_all
    A_bar = A_bar[node_degree_select_indices, :]
    A_bar = A_bar[:, node_degree_select_indices]
    
    
    A = Adj_all
    A_hat = A_pred
    
    #print("Adj_all", Adj_all)
    #print("A_bar", A_bar)
    
    mask = mask_activation(A_bar, H_list, W_list, num_layer)  # last layer (l-1) never used
    
    # =============================================================================
    # def partial_L_partial_A_bar(A, A_hat, A_bar, H, W, mask):   ### equation (1)
    # =============================================================================
    
    n = A_bar.shape[0]
    
    l = num_layer
    
    I_n = sparse_identity(len(A_bar))
    
    U_bar_nn = PermutationRelatedMatrix(n, n)
    
    partial = torch.zeros((A_bar.shape)).to_sparse() ### zero tensor
    
    import time
    
    
    # Start time
    start_time = time.time()
    
    
    #test_interation = few_index_test
    test_interation = len(node_degree_select_indices)
    
    with open("time_collapse.txt", "w") as file:
        for i in range(test_interation):
    
            print('i: ', i)
            for j in range(test_interation):
                print('j: ', j)
                selector_ei = torch.zeros(A_bar.shape[0])
                selector_ei[i] = 1
                selector_ei = selector_ei.unsqueeze(0).to_sparse().T
    
                selector_ej = torch.zeros(A_bar.shape[0])
                selector_ej[j] = 1
                selector_ej = selector_ej.unsqueeze(0).to_sparse()
                selector_ej = selector_ej.T
    
                a_ij = A[i, j]
                hat_a_ij = A_hat[i, j]
                A_hat_ij = hat_a_ij
    
                #%% equation (2)
                # Compute ∂(a_ij log(a_hat_ij))/∂A_bar using equation 2 from your message
    
                partial_aij_A_bar = torch.sparse.mm(torch.sparse.mm( kronecker( I_n, selector_ei.T ) , U_bar_nn)  , ( kronecker( I_n, selector_ej ) ))
    
                # term11 = partial_aij_A_bar @ kronecker( I_n , torch.log(hat_a_ij).unsqueeze(0).unsqueeze(0).to_sparse())
                hat_a_ij = 0
    
                #%% prevent "hat_a_ij = -Inf"
                if hat_a_ij == 0:
                    hat_a_ij += torch.tensor(1e-10)
    
                term11 = torch.sparse.mm(partial_aij_A_bar, ( I_n * torch.log(hat_a_ij) ) )
    
    
                H_li = H_list[-1][i,:]
                H_lj = H_list[-1][j,:]
    
                term12 = ( I_n * a_ij) @  partial_log_term_partial_A_bar(H_li, H_lj, A_bar, A_hat_ij, mask, i, j).detach()
    
                term1 = term11 + term12
    
    
    
    
                #%% partial a_hat_ij to A
                # Compute ∂((1 - a_ij) log(1 - a_hat_ij))/∂A_bar directly
    
                term21 = - partial_aij_A_bar @ ( I_n * torch.log(1 - hat_a_ij)   )
    
                # H_lj_T_to_A_bar
    
                partial_Hi_Hj_T_A_bar_part_1 = H_li_to_A_bar(A_bar, H_list, W_list, I_n, mask, i) @ kronecker( I_n, H_list[l][j,:].unsqueeze(0).T.to_sparse() )
    
    
                partial_Hi_Hj_T_A_bar_part_2 = kronecker( I_n, H_list[l][i,:].unsqueeze(0).to_sparse() ) @ H_lj_T_to_A_bar(A_bar, H_list, W_list, I_n, mask, j)
    
                partial_Hi_Hj_T_A_bar = partial_Hi_Hj_T_A_bar_part_1 + partial_Hi_Hj_T_A_bar_part_2
    
    
                term22 = (I_n* (1 - a_ij)) @  (((-1) * hat_a_ij) * partial_Hi_Hj_T_A_bar)
    
                term2 = term21 + term22
    
                partial -= (term1 + term2)
    
    
                # End time
                end_time = time.time()
    
                # Calculate elapsed time
                elapsed_time = end_time - start_time
    
    #%% Save sensitivity result
# =============================================================================
#     file_path = 'covid_sensitivity_'+str(10*ind)+'_'+str(10*(ind+1))+'.pth'
#     torch.save(partial.to_dense(), file_path)
# =============================================================================





partial_result = partial.to_dense()
partial_result.shape



#%% Load sensitivity result

folder_path = "C://Users//yrt05//Downloads//COVID_pred//COVID_pred"
folder_path = "Data/"

iters = 10
indices_in_A_bar = torch.ones(1)
for ind in range(0,iters):
    minimal_node_indices = selected_indices[min_indices[10*ind: 10*(ind+1)]]
    print(len(minimal_node_indices))
    # minimal_node_indices = selected_indices[min_indices[10:20]]
# =============================================================================
#     ind = 1
# =============================================================================
    
    node_degree_select_indices = minimal_node_indices
    
    
    file_name = 'covid_sensitivity_'+str(10*ind)+'_'+str(10*(ind+1))+'.pth'
    file_path = os.path.join(folder_path, file_name)
    
    # partial_result = torch.load("C://Users//yrt05//Downloads//COVID_pred//COVID_pred//covid_sensitivity.pth")
    partial_result = torch.load(file_path)

    #%% Retrain model
    partial_result_sum = torch.sum(partial_result, dim = 0)
    
    partial_sorted = torch.sort(abs(partial_result_sum), descending = False) ### remove small sensitivity
    
    ### node indices
    # node_indices = node_degree_select[1]

### pick threshold (nodes with highest/lowest degrees)

    pick_threshold = 0.5 # pick first 50% (highest) 

    # len(partial_sorted[1])
    
    node_degree_select_indices_remove = partial_sorted[1][0: int( len(partial_sorted[1]) *pick_threshold) ]
   
    # node_degree_select_indices_remove = torch.sort(node_degree_select_indices_remove, descending = False)[1]
    
    len(node_degree_select_indices_remove)
    
    remove_indices = node_degree_select_indices_remove

    # Assuming A_sub is a subgraph of A_bar and remove_indices are the indices of nodes to be removed from A_sub
    # node_degree_select_indices represent the indices of selected nodes in A_bar
    
    # 1. Create a mask to identify nodes to be removed from A_bar
    mask_remove_from_A_bar = ~torch.isin(torch.arange(len(node_degree_select_indices)), remove_indices)
    
    # 2. Select the corresponding indices of nodes in A_bar
    record_indices = node_degree_select_indices[mask_remove_from_A_bar]

    indices_in_A_bar= torch.cat((indices_in_A_bar, record_indices), dim=0)

indices_in_A_bar = indices_in_A_bar[1:]
print(len(indices_in_A_bar))

print('Unique indices length :', len(torch.unique(indices_in_A_bar)))

#%% Save 'indices_in_A_bar'

# =============================================================================
# file_path = 'covid_indices_in_A_bar.pth'
# torch.save(indices_in_A_bar , file_path)
# =============================================================================

file_path = 'covid_indices_in_A_bar_50_nodes.pth'
torch.save(indices_in_A_bar , file_path)



#
# import networkx as nx
import numpy as np
import torch
import scipy.sparse as sp
from sklearn.neighbors import kneighbors_graph, radius_neighbors_graph
from sklearn.decomposition import PCA  # sklearn PCA is used because PCA in scanpy is not stable.

# from .clustering_func import  mclust_R, louvain
import ot


import pandas as pd
import sklearn.neighbors
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn import metrics

import scanpy as sc

from torch_geometric.data import Data


def refine_label(adata, radius=15, key='label'):
    n_neigh = radius
    new_type = []
    old_type = adata.obs[key].values
    
    #calculate distance
    position = adata.obsm['spatial']
    distance = ot.dist(position, position, metric='euclidean')
    print(distance)
           
    n_cell = distance.shape[0]
    
    for i in range(n_cell):
        vec  = distance[i, :]
        index = vec.argsort()
        neigh_type = []
        for j in range(1, n_neigh+1):
            neigh_type.append(old_type[index[j]])
        max_type = max(neigh_type, key=neigh_type.count)
        new_type.append(max_type)
        
    new_type = [str(i) for i in list(new_type)]    
    adata.obs['label_refined'] = np.array(new_type)
    
    return new_type

#根据细胞类型构建邻接矩阵，没有进行筛选
def generate_celltype_adj_mat(adata, use_rep, include_self=False):

    # sc.pp.neighbors(adata)
    # sc.tl.louvain(adata, resolution=0.8, key_added='expression_louvain_label')

    # new_type = refine_label(adata, key='expression_louvain_label')

    # labels = adata.obs['label_refined']
    # labels = adata.obs['expression_louvain_label']


    adj_mat = np.zeros((len(adata), len(adata)))

    labels = adata.obs[use_rep]
    
    # labels = adata.obs['celltype']

    for i in range(len(adata)):
        n_neighbors = np.where(labels == labels[i])
        adj_mat[i, n_neighbors] = 1

    if not include_self:
        x, y = np.diag_indices_from(adj_mat)
        adj_mat[x, y] = 0

    adj_mat = adj_mat + adj_mat.T
    adj_mat = adj_mat > 0
    adj_mat = adj_mat.astype(np.int64)
    
    return adj_mat


def generate_coordinate_adj_mat_by_kneighbors(adata, include_self=False, n=6):

    assert 'spatial' in adata.obsm, 'AnnData object should provided spatial information'

    # dist = metrics.pairwise_distances(adata.obsm['spatial'])
    coor = pd.DataFrame(adata.obsm['spatial'])
    coor.index = adata.obs.index
    coor.columns = ['imagerow', 'imagecol']

    nbrs = sklearn.neighbors.NearestNeighbors(n_neighbors=n+1).fit(coor)
    distances, indices = nbrs.kneighbors(coor)

    adata.uns['Spatial_Net'] = indices

    adj_mat = np.zeros((len(adata), len(adata)))
    for i in range(len(adata)):
        n_neighbors = indices[i]
        # n_neighbors = np.argsort(dist[i, :])[:n+1]
        adj_mat[i, n_neighbors] = 1

    if not include_self:
        x, y = np.diag_indices_from(adj_mat)
        adj_mat[x, y] = 0

    adj_mat = adj_mat + adj_mat.T
    adj_mat = adj_mat > 0
    adj_mat = adj_mat.astype(np.int64)

    return adj_mat


def generate_coordinate_adj_mat_by_radius(adata, include_self=False, n_radius=150):
    assert 'spatial' in adata.obsm, 'AnnData object should provided spatial information'

    coor = pd.DataFrame(adata.obsm['spatial'])
    coor.index = adata.obs.index
    coor.columns = ['imagerow', 'imagecol']

    # nbrs = sklearn.neighbors.NearestNeighbors(n_neighbors=15+1).fit(coor)
    # distances, indices = nbrs.kneighbors(coor)
    # KNN_list = []
    # for it in range(indices.shape[0]):
    #     KNN_list.append(pd.DataFrame(zip([it]*indices.shape[1],indices[it,:], distances[it,:])))


    nbrs = sklearn.neighbors.NearestNeighbors(radius=n_radius).fit(coor)
    distances, indices = nbrs.radius_neighbors(coor, return_distance=True)

    adata.uns['Spatial_Net'] = indices


    adj_mat = np.zeros((len(adata), len(adata)))

    for i in range(len(adata)):
        adj_mat[i, indices[i]] = 1

    if not include_self:
        x, y = np.diag_indices_from(adj_mat)
        adj_mat[x, y] = 0

    adj_mat = adj_mat + adj_mat.T
    adj_mat = adj_mat > 0
    adj_mat = adj_mat.astype(np.int64)
    
    return adj_mat



def generate_celltype_coordinate_adj_mat(adata, include_self=False, n=6):
    import sklearn.neighbors

    from sklearn import metrics
    assert 'spatial' in adata.obsm, 'AnnData object should provided spatial information'

    # dist = metrics.pairwise_distances(adata.obsm['spatial'])

    # coor = pd.DataFrame(adata.obsm['spatial'])
    # coor.index = adata.obs.index
    # coor.columns = ['imagerow', 'imagecol']


    # nbrs = sklearn.neighbors.NearestNeighbors(n_neighbors=15+1).fit(coor)
    # distances, indices = nbrs.kneighbors(coor)
    # KNN_list = []
    # for it in range(indices.shape[0]):
    #     KNN_list.append(pd.DataFrame(zip([it]*indices.shape[1],indices[it,:], distances[it,:])))

    # nbrs = sklearn.neighbors.NearestNeighbors(n_neighbors=15+1).fit(coor)
    # distances, indices = nbrs.kneighbors(coor)

    indices = adata.uns['Spatial_Net']


    # nbrs = sklearn.neighbors.NearestNeighbors(radius=200).fit(coor)
    # distances, indices = nbrs.radius_neighbors(coor, return_distance=True)
    # print(distances[0])
    # print(indices[0])

    # sorted_indices = np.argsort(distances[0])
    # sorted_distances = distances[0][sorted_indices]
    # sorted_neighbors = indices[0][sorted_indices]
    # print(sorted_distances)
    # print(sorted_neighbors)


    adj_mat = np.zeros((len(adata), len(adata)))

    # sc.pp.neighbors(adata)
    # sc.tl.louvain(adata, resolution=0.8, key_added='expression_louvain_label')


    # labels = adata.obs['expression_louvain_label']

    labels = adata.obs['layer_guess']

    for i in range(len(adata)):
        n_neighbors = indices[i]
        for k in range(len(n_neighbors)):
            if(labels[i] == labels[n_neighbors[k]]):
                adj_mat[i, n_neighbors[k]] = 1


    if not include_self:
        x, y = np.diag_indices_from(adj_mat)
        adj_mat[x, y] = 0

    adj_mat = adj_mat + adj_mat.T
    adj_mat = adj_mat > 0
    adj_mat = adj_mat.astype(np.int64)
    
    return adj_mat



def generate_feature_adj_mat(adata, use_pre, k=20, mode= "connectivity", metric="correlation", include_self=False):
        
    feature_graph = kneighbors_graph(adata.obsm[use_pre], k, mode=mode, metric=metric, include_self=include_self)
    # feature_graph = radius_neighbors_graph(adata.obsm[use_pre], radius=1, mode=mode, metric=metric, include_self=include_self)

    
    adj_mat = feature_graph.toarray()

    adj_mat = adj_mat + adj_mat.T
    adj_mat = adj_mat > 0
    adj_mat = adj_mat.astype(np.int64)
    
    return adj_mat






def generate_celltype_adj_mat1(adata, include_self=False, n=6):
    from sklearn import metrics
    assert 'spatial' in adata.obsm, 'AnnData object should provided spatial information'

    # dist = metrics.pairwise_distances(adata.obsm['spatial'])

    coor = pd.DataFrame(adata.obsm['spatial'])
    coor.index = adata.obs.index
    coor.columns = ['imagerow', 'imagecol']


    # nbrs = sklearn.neighbors.NearestNeighbors(n_neighbors=15+1).fit(coor)
    # distances, indices = nbrs.kneighbors(coor)
    # KNN_list = []
    # for it in range(indices.shape[0]):
    #     KNN_list.append(pd.DataFrame(zip([it]*indices.shape[1],indices[it,:], distances[it,:])))


    nbrs = sklearn.neighbors.NearestNeighbors(radius=300).fit(coor)
    distances, indices = nbrs.radius_neighbors(coor, return_distance=True)


    adj_mat = np.zeros((len(adata), len(adata)))

    labels = adata.obs['layer_guess']

    for i in range(len(adata)):
        n_neighbors = []

        sorted_index = np.argsort(-distances[i])
        sorted_distances = distances[i][sorted_index]
        sorted_neighbors = indices[i][sorted_index]

        for k in range(len(sorted_neighbors)):
            if(labels[i] == labels[sorted_neighbors[k]]):
                n_neighbors.append(sorted_neighbors[k])
        adj_mat[i,n_neighbors[:6]] = 1
        adj_mat[i,n_neighbors[5:]] = 1

    if not include_self:
        x, y = np.diag_indices_from(adj_mat)
        adj_mat[x, y] = 0

    adj_mat = adj_mat + adj_mat.T
    adj_mat = adj_mat > 0
    adj_mat = adj_mat.astype(np.int64)
    
    return adj_mat

##### generate n
def generate_adj_mat(adata, include_self=False, n=6):
    from sklearn import metrics
    assert 'spatial' in adata.obsm, 'AnnData object should provided spatial information'

    dist = metrics.pairwise_distances(adata.obsm['spatial'])

    # sample_name = list(adata.uns['spatial'].keys())[0]
    # scalefactors = adata.uns['spatial'][sample_name]['scalefactors']
    # adj_mat = dist <= scalefactors['fiducial_diameter_fullres'] * (n+0.2)
    # adj_mat = adj_mat.astype(int)

    # n_neighbors = np.argpartition(dist, n+1, axis=1)[:, :(n+1)]
    # adj_mat = np.zeros((len(adata), len(adata)))
    # for i in range(len(adata)):
    #     adj_mat[i, n_neighbors[i, :]] = 1

    adj_mat = np.zeros((len(adata), len(adata)))
    for i in range(len(adata)):
        n_neighbors = np.argsort(dist[i, :])[:n+1]
        adj_mat[i, n_neighbors] = 1

    if not include_self:
        x, y = np.diag_indices_from(adj_mat)
        adj_mat[x, y] = 0

    adj_mat = adj_mat + adj_mat.T
    adj_mat = adj_mat > 0
    adj_mat = adj_mat.astype(np.int64)
    
    return adj_mat


def generate_adj_mat_radius(adata, include_self=False, n=6):
    
    coor = pd.DataFrame(adata.obsm['spatial'])
    coor.index = adata.obs.index
    coor.columns = ['imagerow', 'imagecol']


    # nbrs = sklearn.neighbors.NearestNeighbors(n_neighbors=15+1).fit(coor)
    # distances, indices = nbrs.kneighbors(coor)
    # KNN_list = []
    # for it in range(indices.shape[0]):
    #     KNN_list.append(pd.DataFrame(zip([it]*indices.shape[1],indices[it,:], distances[it,:])))


    nbrs = sklearn.neighbors.NearestNeighbors(radius=200).fit(coor)
    distances, indices = nbrs.radius_neighbors(coor, return_distance=True)
    KNN_list = []
    for it in range(indices.shape[0]):
        KNN_list.append(pd.DataFrame(zip([it]*indices[it].shape[0], indices[it], distances[it])))

    KNN_df = pd.concat(KNN_list)
    KNN_df.columns = ['Cell1', 'Cell2', 'Distance']

    Spatial_Net = KNN_df.copy()
    Spatial_Net = Spatial_Net.loc[Spatial_Net['Distance']>0,]
    id_cell_trans = dict(zip(range(coor.shape[0]), np.array(coor.index), ))
    Spatial_Net['Cell1'] = Spatial_Net['Cell1'].map(id_cell_trans)
    Spatial_Net['Cell2'] = Spatial_Net['Cell2'].map(id_cell_trans)
    
    adata.uns['Spatial_Net'] = Spatial_Net

    G_df = adata.uns['Spatial_Net'].copy()
    cells = np.array(adata.obs_names)
    cells_id_tran = dict(zip(cells, range(cells.shape[0])))
    G_df['Cell1'] = G_df['Cell1'].map(cells_id_tran)
    G_df['Cell2'] = G_df['Cell2'].map(cells_id_tran)

    G = sp.coo_matrix((np.ones(G_df.shape[0]), (G_df['Cell1'], G_df['Cell2'])), shape=(adata.n_obs, adata.n_obs))

    print("**********Radius Graph************")

    return G

def generate_adj_mat_celtype(adata, include_self=False, n=6):
    
    coor = pd.DataFrame(adata.obsm['spatial'])
    coor.index = adata.obs.index
    coor.columns = ['imagerow', 'imagecol']

    nbrs = sklearn.neighbors.NearestNeighbors(radius=200).fit(coor)
    distances, indices = nbrs.radius_neighbors(coor, return_distance=True)
    KNN_list = []
    for it in range(indices.shape[0]):
        KNN_list.append(pd.DataFrame(zip([it]*indices[it].shape[0], indices[it], distances[it])))

    KNN_df = pd.concat(KNN_list)
    KNN_df.columns = ['Cell1', 'Cell2', 'Distance']

    Spatial_Net = KNN_df.copy()
    Spatial_Net = Spatial_Net.loc[Spatial_Net['Distance']>0,]
    id_cell_trans = dict(zip(range(coor.shape[0]), np.array(coor.index), ))
    Spatial_Net['Cell1'] = Spatial_Net['Cell1'].map(id_cell_trans)
    Spatial_Net['Cell2'] = Spatial_Net['Cell2'].map(id_cell_trans)
    
    adata.uns['Spatial_Net'] = Spatial_Net

    Graph_df = adata.uns['Spatial_Net'].copy()

    # sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata)
    sc.tl.louvain(adata, resolution=0.2, key_added='expression_louvain_label')
    # sc.tl.leiden(adata, resolution=0.2, key_added='expression_louvain_label')


    # louvain(adata, n_clusters=7, use_rep='X_pca', key_added='expression_louvain_label')

    label = adata.obs['expression_louvain_label']

    print('------Pruning the graph...')
    print('%d edges before pruning.' %Graph_df.shape[0])

    pro_labels_dict = dict(zip(list(label.index), label))
    Graph_df['Cell1_label'] = Graph_df['Cell1'].map(pro_labels_dict)
    Graph_df['Cell2_label'] = Graph_df['Cell2'].map(pro_labels_dict)
    Graph_df = Graph_df.loc[Graph_df['Cell1_label']==Graph_df['Cell2_label'],]

    print('%d edges after pruning.' %Graph_df.shape[0])

    adata_Vars =  adata[:, adata.var['highly_variable']]
    X = pd.DataFrame(adata_Vars.X.toarray()[:, ], index=adata_Vars.obs.index, columns=adata_Vars.var.index)
    cells = np.array(X.index)
    cells_id_tran = dict(zip(cells, range(cells.shape[0])))

    prune_G_df = Graph_df
    prune_G_df['Cell1'] = prune_G_df['Cell1'].map(cells_id_tran)
    prune_G_df['Cell2'] = prune_G_df['Cell2'].map(cells_id_tran)

    prune_G = sp.coo_matrix((np.ones(prune_G_df.shape[0]), (prune_G_df['Cell1'], prune_G_df['Cell2'])), shape=(adata.n_obs, adata.n_obs))

    return prune_G

def generate_adj_mat_celtype1(adata, include_self=False, n=6):
    
    coor = pd.DataFrame(adata.obsm['spatial'])
    coor.index = adata.obs.index
    coor.columns = ['imagerow', 'imagecol']


    # nbrs = sklearn.neighbors.NearestNeighbors(n_neighbors=30+1).fit(coor)
    # distances, indices = nbrs.kneighbors(coor)
    # KNN_list = []
    # for it in range(indices.shape[0]):
    #     KNN_list.append(pd.DataFrame(zip([it]*indices.shape[1],indices[it,:], distances[it,:])))


    nbrs = sklearn.neighbors.NearestNeighbors(radius=150).fit(coor)
    distances, indices = nbrs.radius_neighbors(coor, return_distance=True)
    KNN_list = []
    for it in range(indices.shape[0]):
        KNN_list.append(pd.DataFrame(zip([it]*indices[it].shape[0], indices[it], distances[it])))

    KNN_df = pd.concat(KNN_list)
    KNN_df.columns = ['Cell1', 'Cell2', 'Distance']

    Spatial_Net = KNN_df.copy()
    Spatial_Net = Spatial_Net.loc[Spatial_Net['Distance']>0,]
    id_cell_trans = dict(zip(range(coor.shape[0]), np.array(coor.index), ))
    Spatial_Net['Cell1'] = Spatial_Net['Cell1'].map(id_cell_trans)
    Spatial_Net['Cell2'] = Spatial_Net['Cell2'].map(id_cell_trans)
    
    adata.uns['Spatial_Net'] = Spatial_Net

    Graph_df = adata.uns['Spatial_Net'].copy()

    # ******************label***********************
    category = adata.obs['layer_guess'].unique()
    category_id_tran = dict(zip(category, range(category.shape[0])))

    cell_type = pd.DataFrame(adata.obs['layer_guess'])
    cell_type.index = adata.obs.index
    cell_type['layer_guess'] = cell_type['layer_guess'].map(category_id_tran)

    adata.obs['pre_labela'] = cell_type
    label = adata.obs['pre_labela']

    print('------Pruning the graph...')
    print('%d edges before pruning.' %Graph_df.shape[0])

    pro_labels_dict = dict(zip(list(label.index), label))
    Graph_df['Cell1_label'] = Graph_df['Cell1'].map(pro_labels_dict)
    Graph_df['Cell2_label'] = Graph_df['Cell2'].map(pro_labels_dict)
    Graph_df = Graph_df.loc[Graph_df['Cell1_label']==Graph_df['Cell2_label'],]

    print('%d edges after pruning.' %Graph_df.shape[0])

    adata_Vars =  adata[:, adata.var['highly_variable']]
    X = pd.DataFrame(adata_Vars.X.toarray()[:, ], index=adata_Vars.obs.index, columns=adata_Vars.var.index)
    cells = np.array(X.index)
    cells_id_tran = dict(zip(cells, range(cells.shape[0])))

    prune_G_df = Graph_df
    prune_G_df['Cell1'] = prune_G_df['Cell1'].map(cells_id_tran)
    prune_G_df['Cell2'] = prune_G_df['Cell2'].map(cells_id_tran)

    prune_G = sp.coo_matrix((np.ones(prune_G_df.shape[0]), (prune_G_df['Cell1'], prune_G_df['Cell2'])), shape=(adata.n_obs, adata.n_obs))

    return prune_G




def generate_adj_mat_1(adata, max_dist):
    from sklearn import metrics
    assert 'spatial' in adata.obsm, 'AnnData object should provided spatial information'

    dist = metrics.pairwise_distances(adata.obsm['spatial'], metric='euclidean')
    adj_mat = dist < max_dist
    adj_mat = adj_mat.astype(np.int64)
    return adj_mat

##### normalze graph
def sparse_mx_to_torch_sparse_tensor(sparse_mx):
    """Convert a scipy sparse matrix to a torch sparse tensor."""
    sparse_mx = sparse_mx.tocoo().astype(np.float32)
    indices = torch.from_numpy(np.vstack((sparse_mx.row, sparse_mx.col)).astype(np.int64))
    values = torch.from_numpy(sparse_mx.data)
    shape = torch.Size(sparse_mx.shape)
    return torch.sparse.FloatTensor(indices, values, shape)

def preprocess_graph(adj):
    adj_ = adj + sp.eye(adj.shape[0])
    rowsum = np.array(adj_.sum(1))
    degree_mat_inv_sqrt = sp.diags(np.power(rowsum, -0.5).flatten())
    adj_normalized = adj_.dot(degree_mat_inv_sqrt).transpose().dot(degree_mat_inv_sqrt).tocoo()
    return sparse_mx_to_torch_sparse_tensor(adj_normalized)


def mask_generator(adj_label, N=1):
    idx = adj_label.indices()
    cell_num = adj_label.size()[0]

    list_non_neighbor = []
    for i in range(0, cell_num):
        neighbor = idx[1, torch.where(idx[0, :] == i)[0]]
        n_selected = len(neighbor) * N

        # non neighbors
        total_idx = torch.range(0, cell_num-1, dtype=torch.float32)
        non_neighbor = total_idx[~torch.isin(total_idx, neighbor)]
        indices = torch.randperm(len(non_neighbor), dtype=torch.float32)
        random_non_neighbor = indices[:n_selected]
        list_non_neighbor.append(random_non_neighbor)

    x = adj_label.indices()[0]
    y = torch.concat(list_non_neighbor)

    indices = torch.stack([x, y])
    indices = torch.concat([adj_label.indices(), indices], axis=1)

    value = torch.concat([adj_label.values(), torch.zeros(len(x), dtype=torch.float32)])
    adj_mask = torch.sparse_coo_tensor(indices, value)

    return adj_mask


def graph_computing(pos, n):
    from scipy.spatial import distance
    list_x = []
    list_y = []
    list_value = []

    for node_idx in range(len(pos)):
        tmp = pos[node_idx, :].reshape(1, -1)
        distMat = distance.cdist(tmp, pos, 'euclidean')
        res = distMat.argsort()
        # tmpdist = distMat[0, res[0][1:params.k + 1]]
        for j in np.arange(1, n + 1):
            list_x += [node_idx, res[0][j]]
            list_y += [res[0][j], node_idx]
            list_value += [1, 1]

    adj = sp.csr_matrix((list_value, (list_x, list_y)))
    adj = adj >= 1
    adj = adj.astype(np.float32)
    return adj


def graph_construction(adata, n=6, n_radius=150, feat_use='X_pca', celltype_use='layer_guess',k=20, dmax=50, mode='KNN'):
    if mode == 'KNN':
        adj_m1 = generate_coordinate_adj_mat_by_kneighbors(adata, include_self=False, n=n)
        # adj_m1 = generate_celltype_coordinate_adj_mat(adata)
        #3639 x 3639
        # adj_m1 = graph_computing(adata.obsm['spatial'], n=n)
    elif mode == 'Radius':
        adj_m1 = generate_coordinate_adj_mat_by_radius(adata, include_self=False, n_radius=n_radius)
    elif mode == 'Feature':
        adj_m1 = generate_feature_adj_mat(adata, use_pre=feat_use, k=k, include_self=False)
    elif mode == 'Celltype':
        adj_m1 = generate_celltype_adj_mat(adata, use_rep=celltype_use)
        # adj_m1 = generate_celltype_coordinate_adj_mat(adata)
    else:
        adj_m1 = generate_adj_mat_1(adata, dmax)
    adj_m1 = sp.coo_matrix(adj_m1)
    # adj_m1 = generate_adj_mat(adata, include_self=False, n=n)

    # Store original adjacency matrix (without diagonal entries) for later
    adj_m1 = adj_m1 - sp.dia_matrix((adj_m1.diagonal()[np.newaxis, :], [0]), shape=adj_m1.shape)
    adj_m1.eliminate_zeros()

    # Some preprocessing
    adj_norm_m1 = preprocess_graph(adj_m1)
    adj_m1 = adj_m1 + sp.eye(adj_m1.shape[0])
    # adj_label_m1 = torch.FloatTensor(adj_label_m1.toarray())

    adj_m1 = adj_m1.tocoo()
    shape = adj_m1.shape
    values = adj_m1.data
    indices = np.stack([adj_m1.row, adj_m1.col])
    adj_label_m1 = torch.sparse_coo_tensor(indices, values, shape)


    if(float((adj_m1.shape[0] * adj_m1.shape[0] - (adj_m1.sum())) * 2) == 0):
        norm_m1 = adj_m1.shape[0] * adj_m1.shape[0] / float((adj_m1.shape[0] * adj_m1.shape[0] - (adj_m1.sum()-0.0001)) * 2)
    else:
        norm_m1 = adj_m1.shape[0] * adj_m1.shape[0] / float((adj_m1.shape[0] * adj_m1.shape[0] - (adj_m1.sum())) * 2)


    # # generate random mask
    # adj_mask = mask_generator(adj_label_m1.to_sparse(), N)

    graph_dict = {
        "adj_norm": adj_norm_m1,
        "adj_label": adj_label_m1.coalesce(),
        "norm_value": norm_m1,
        # "mask": adj_mask
    }

    return graph_dict


def graph_construction_celltype(adata, n=6, dmax=50, mode='KNN'):
    # if mode == 'KNN':
    #     adj_m1 = generate_adj_mat(adata, include_self=False, n=n)
    #     #3639 x 3639
    #     # adj_m1 = graph_computing(adata.obsm['spatial'], n=n)
    # else:
    #     adj_m1 = generate_adj_mat_1(adata, dmax)
    # adj_m1 = sp.coo_matrix(adj_m1)

    adj_m1 = generate_adj_mat_celtype1(adata, include_self=False, n=n)

    # Store original adjacency matrix (without diagonal entries) for later
    adj_m1 = adj_m1 - sp.dia_matrix((adj_m1.diagonal()[np.newaxis, :], [0]), shape=adj_m1.shape)
    adj_m1.eliminate_zeros()

    # Some preprocessing
    adj_norm_m1 = preprocess_graph(adj_m1)
    adj_m1 = adj_m1 + sp.eye(adj_m1.shape[0])
    # adj_label_m1 = torch.FloatTensor(adj_label_m1.toarray())

    adj_m1 = adj_m1.tocoo()
    shape = adj_m1.shape
    values = adj_m1.data
    indices = np.stack([adj_m1.row, adj_m1.col])
    adj_label_m1 = torch.sparse_coo_tensor(indices, values, shape)

    norm_m1 = adj_m1.shape[0] * adj_m1.shape[0] / float((adj_m1.shape[0] * adj_m1.shape[0] - adj_m1.sum()) * 2)

    # # generate random mask
    # adj_mask = mask_generator(adj_label_m1.to_sparse(), N)

    graph_dict = {
        "adj_norm": adj_norm_m1,
        "adj_label": adj_label_m1.coalesce(),
        "norm_value": norm_m1,
        # "mask": adj_mask
    }

    return graph_dict





def combine_graph_dict(dict_1, dict_2):
    # TODO add adj_org
    tmp_adj_norm = torch.block_diag(dict_1['adj_norm'].to_dense(), dict_2['adj_norm'].to_dense())
    tmp_adj_label = torch.block_diag(dict_1['adj_label'].to_dense(), dict_2['adj_label'].to_dense())
    graph_dict = {
        "adj_norm": tmp_adj_norm.to_sparse(),
        "adj_label": tmp_adj_label.to_sparse(),
        "norm_value": np.mean([dict_1['norm_value'], dict_2['norm_value']])
    }
    return graph_dict
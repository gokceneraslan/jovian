# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.2.4
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # jovian - 04 Regulatory Networks

# %%
import scanpy as sc

# Plotting
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

# numpy et al.
import numpy as np
import scipy.sparse as sp
import pandas as pd

from sklearn.model_selection import GridSearchCV

from inverse_covariance import (
    QuicGraphicalLasso,
    QuicGraphicalLassoCV,
    QuicGraphicalLassoEBIC,
    AdaptiveGraphicalLasso,
    ModelAverage,
)

# R integration
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.vectors import StrVector, FloatVector, ListVector, IntVector
from rpy2.robjects import r
import rpy2.robjects as ro
import anndata2ri # scipy.sparse + AnnData support

numpy2ri.activate()
pandas2ri.activate()
anndata2ri.activate()

# %%
# %load_ext rpy2.ipython

# %% [markdown]
# ## Parameters

# %%
# %store -r species
# %store -r save_filename
# %store -r save_filename_group

# %%
metric = 'log_likelihood'
num_folds = 3
n_jobs = 3

# %% [markdown]
# ## Read dataset

# %%
adata = sc.read('adata-brain.h5ad')
adata

# %%
sc.pp.highly_variable_genes(adata, n_top_genes=200, subset=True)
adata

# %%
sc.pl.highly_variable_genes(adata)

# %%
adata = adata[:2000].copy()
adata = adata[:, (adata.X!=0).sum(0).A1 > 20].copy()
adata


# %% [markdown]
# ## mgm

# %%
def generate_mgm_input(adata, keys=None, gene_poisson=False):
    adata._sanitize()
    if keys is None: keys = []

    mat = []
    nlevels = np.ones((adata.n_vars,)).tolist()
    gene_type = 'p' if gene_poisson else 'g'
    types = np.repeat(gene_type, adata.n_vars).tolist()

    if keys:
        obs = adata.obs[keys]

        for key in keys:
            if obs[key].dtype.name == 'category':
                mat.append((obs[key].cat.codes+1).values[:, None])
                nlevels.append(len(obs[key].cat.categories))
                types.append('c')
            else:
                mat.append(obs[key].values[:, None])
                nlevels.append(1)
                types.append('g')
    if mat:
        mat = np.concatenate(mat, axis=1)
        if sp.issparse(adata.X):
            mat = sp.hstack([adata.X, sp.coo_matrix(mat)])
        else:
            mat = np.hstack([adata.X, mat])
    else:
        mat = adata.X
    labels = adata.var_names.tolist() + keys
    node_colors = sns.color_palette('tab10', n_colors=len(keys)).as_hex()
    
    return mat, nlevels, types, labels, (adata.n_vars * ['#ffffff']) + node_colors


def fit_mgm(adata, keys=None, gene_poisson=False, lambdaGam=0.5, **kwds):

    mat, nlevels, types, labels, node_colors = generate_mgm_input(adata, keys=keys, gene_poisson=gene_poisson)

    mgm = importr('mgm')
    fit = mgm.mgm(data=r['as.matrix'](mat), 
                  type=StrVector(types), 
                  level=IntVector(nlevels), 
                  saveModels=False,
                  lambdaGam=lambdaGam,
                  **kwds)
    
    ret = {}

    ret['fit'] = fit
    ret['wadj'] = fit.rx2('pairwise').rx2('wadj')
    ret['signs'] = fit.rx2('pairwise').rx2('signs')
    ret['edge_colors'] = fit.rx2('pairwise').rx2('edgecolor')
    ret['node_colors'] = node_colors
    ret['node_labels'] = labels
    
    return ret



# %%
mgm = fit_mgm(adata, ['predicted_cell_types'])

# %%
mgm['signs'].flatten()[np.isfinite(mgm['signs'].flatten())]

# %%
import networkx as nx
import hvplot.networkx as hvnx
import igraph as ig
from scanpy._utils import get_igraph_from_adjacency


def plot_rene(fit, remove_key_nodes=True, layout='fr', width=1000, height=1000):

    adj = fit['wadj']
    signs = fit['signs']
    node_labels = fit['node_labels']
    node_colors = fit['node_colors']
    edge_colors = fit['edge_colors']
    
    n_keys = len(set(node_colors)) - 1
    
    if remove_key_nodes:
        node_labels = node_labels[:-n_keys]
        node_colors = ['lightgray' if x else 'white' for x in (adj[-n_keys:, :-n_keys].sum(0) != 0)]
        adj = adj[:-n_keys, :-n_keys]
        signs = signs[:-n_keys, :-n_keys]

    g = get_igraph_from_adjacency(adj, directed=False)
    l = {i: v for i, v in enumerate(node_labels)}

    G = nx.from_numpy_matrix(adj*10)
    
    edge_colors = [(signs*adj)[x, y] for x,y in G.edges] 
    
    G = nx.relabel_nodes(G, l)
    layout = g.layout(layout)
    pos = {i: v for i, v in zip(node_labels, layout)}

    return hvnx.draw(G,
              pos=pos,
              font_size='7pt', 
              node_size=800,
              edge_width='weight',
              #edge_cmap='coolwarm',
              #edge_color=edge_colors,
              node_color=node_colors,
              with_labels=True, 
              width=width,
              height=height)    
    
plot_rene(mgm, remove_key_nodes=True)

# %%
fit = mgm['fit']

# %% {"magic_args": "-i labels -i fit -u in -w 10 -h 10 -r 300", "language": "R"}
#
# FactorGraph(fit, labels=unlist(labels), layout='sugiyama', PairwiseAsEdge=T, repulsion=1.3, shapeSizes=c(1, 0.05))

# %% {"magic_args": "-i labels -u in -w 20 -h 20 -r 200", "language": "R"}
#
# library(qgraph)
#
# ri = rowSums(fit$pairwise$wadj) > 0.1
# ci = colSums(fit$pairwise$wadj) > 0.1
# labels = unlist(labels)[ri]
#
# qgraph(fit$pairwise$wadj[ri, ci],
#        edge.color = fit$pairwise$edgecolor[ri, ci],
#        layout = "spring",
#        repulsion=1.3,
#        vsize=1,
#        labels = labels)

# %% [markdown]
# ## qgraph

# %%
# %%time
# %%R -i corr -i labels -i nobs -u in -w 20 -h 20 -r 200

library(qgraph)

labels = unlist(labels)
fit = EBICglasso(corr, n=nobs, gamma=0.25, checkPD=F)
qgraph(fit, sampleSize=nobs, layout='spring', labels=labels, repulsion=1.2)
#qgraph(corr, graph='pcor', sampleSize=nobs, layout='spring', labels=labels, repulsion=1.2,  minimum='sig', bonf=T)

# %% [markdown]
# ## gLASSO

# %% {"magic_args": "-i X -o fit", "language": "R"}
#
# library(CVglasso)
# library(glasso)
#
# fit = CVglasso(X, cores=2, trace=T)
#
# plot(fit)

# %%
(wi != 0).sum()

# %% [markdown]
# ## Plot

# %%
from igraph import *

# %%
wi = (wi != 0).astype(float)
wi = wi * wi.T
np.fill_diagonal(wi, 0)

# %%
idx = wi.sum(1) > 0

# %%
genes = adata.var_names[idx]
genes

# %%
wi = wi[np.ix_(idx, idx)]

# %%
g = sc._utils.get_igraph_from_adjacency(wi)

# %%
g.vs['gene'] = genes

# %%
g.layout('fr')

# %%
g.ecount()

# %%

# %%
sc.pl.umap(adata, color=['Krt18', 'Dapl1', 'Defb11'], cmap='Reds')

# %%
g.vcount(), g.ecount()

# %%
plot(g, vertex_size=3, edge_width=0.1, vertex_label=g.vs['gene'], vertex_label_size=4, vertex_label_color='blue')


# %% {"jupyter": {"source_hidden": true}}
def _link_cluster_from_adj(A, resolution=1.0, directed=True, seed=0):
    from scanpy._utils import get_igraph_from_adjacency
    import leidenalg as la
        
    gi = get_igraph_from_adjacency(A, directed=directed)
    print(gi.vcount(), gi.ecount(), flush=True)
    
    gi_line = gi.linegraph()
    print(gi_line.vcount(), gi_line.ecount(), flush=True)
    
    partition = la.find_partition(gi_line, 
                                  la.RBConfigurationVertexPartition, 
                                  resolution_parameter=resolution,
                                  seed=seed)

    n_clusters = len(set(partition.membership))
    print(n_clusters, flush=True)
    
    node2cluster = {}
    for i in trange(gi.ecount()):
        e = gi.es[i]
        cluster_id = partition.membership[i]

        node2cluster.setdefault(e.source, set()).add(cluster_id)
        node2cluster.setdefault(e.target, set()).add(cluster_id)

    return node2cluster, n_clusters

def link_cluster(adata, resolution=1.0, directed=True, bipartite=False, uncertainty=True, seed=0):
    
    if not bipartite:
        A = (adata.uns['neighbors']['connectivities'] != 0).astype(int)
    else:
        ncells, ngenes = adata.n_obs, adata.n_vars
        nnodes = ncells + ngenes
        
        A = sp.lil_matrix((nnodes, nnodes))
        A[:ncells, ncells:] = adata.X
        A = A.tocoo().tocsr()

    node2cluster, n_clusters = _link_cluster_from_adj(A,
                                                      resolution=resolution,
                                                      directed=directed,
                                                      seed=seed)   
    if bipartite:
        n2c_cells = {}
        n2c_genes = {}

        for k, v in node2cluster.items():
            if k < adata.n_obs:
                n2c_cells[k] = v
            else:
                n2c_genes[adata.var_names[k-adata.n_obs]] = v
                #TODO: store n2c_genes
        node2cluster = n2c_cells

    for i in range(n_clusters):
        if f'cluster_{i}' in adata.obs:
            del adata.obs[f'cluster_{i}']
        adata.obs[f'cluster_{i}'] = np.zeros(adata.n_obs, dtype=int)
    adata.obs[f'link_cluster_uncertainty'] = np.zeros(adata.n_obs, dtype=int)

    for k, v in node2cluster.items():
        for cluster_id in v:
            adata.obs[f'cluster_{cluster_id}'].iloc[k] = 1./len(node2cluster[k]) if uncertainty else 1.
        adata.obs[f'link_cluster_uncertainty'].iloc[k] = len(node2cluster[k])
        
    adata.obs[f'link_cluster_uncertainty'] = np.log10(adata.obs[f'link_cluster_uncertainty'])
    
    adata.uns['link_cluster'] = dict(n_clusters=n_clusters,
                                     cluster_names=[f'cluster_{i}' for i in range(n_clusters)])

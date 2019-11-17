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
# # jovian - 04 Velocyto

# %%
import scanpy as sc
import scvelo as scv

# Plotting
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

# numpy et al.
import numpy as np
import scipy.sparse as sp
import pandas as pd

# R integration
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, numpy2ri, r
import rpy2.robjects as ro
import anndata2ri

from pathlib import Path
import math
from tqdm.auto import tqdm
import warnings

# %%
# %load_ext rpy2.ipython

# %%
sc.set_figure_params(dpi=100, fontsize=12)
matplotlib.rcParams['font.sans-serif'] = matplotlib.rcParamsDefault['font.sans-serif']

sc.settings.verbosity = 'hint'

# %%
plt.rcParams['figure.figsize'] = (7, 7)

# %% [markdown]
# ## Parameters

# %% [markdown]
# ## Read processed AnnData

# %%
adata = sc.read('adata-velocyto.h5ad')
adata

# %%
adata.var_names.str.startswith('Mir').sum()

# %% [markdown]
# ## Transfer manual annotations 

# %%
labels = []
for p in ('P1', 'P7', 'P21'):
    df = pd.read_excel('metadata/cell-annotations.xlsx', sheet_name=p)[['Unnamed: 0', 'Consolidated cell type']]
    df['Consolidated cell type'] = [x.strip() for x in df['Consolidated cell type']]
    #df.Remove = [x.strip() for x in df.Remove]
    df.rename(columns={'Unnamed: 0': 'sample_leiden',
                       'Consolidated cell type': 'Cell type'}, inplace=True)
    df.sample_leiden = [x.split()[1] for x in df.sample_leiden]
    labels.append(df.assign(postnatal_age=int(p[1:])))
labels = pd.concat(labels, 0)

labels

# %%
#adata.obs.drop(['Consolidated cell type'], axis=1, inplace=True)

# %%
#del adata.uns['Cell type_colors']

# %%
adata.obs = adata.obs.merge(labels, how='left', on=['sample_leiden', 'postnatal_age']).set_index(adata.obs.index)

# %%
adata.obs

# %%
adata.write('adata-velocyto.h5ad')

# %% [markdown]
# labels = pd.read_csv('metadata/cell-type-annotations.csv', index_col=0)
# new_ind = []
# for l in labels.index:
#     barcode, sample = l.split('-')
#     new_ind.append(f'{sample}:{barcode}x-{sample}')
# labels.index = new_ind
# labels.head()

# %%
#adata.obs = adata.obs.join(labels)

# %%
sc.pl.umap(adata, color='Cell type', size=1)

# %%
sc.pl.dotplot(adata, var_names=adata.var_names[adata.var_names.str.startswith('Mir')], groupby='Cell type', standard_scale='var', mean_only_expressed=True)

# %%
sc.pl.umap(adata, color=['Mir9-3hg', 'Mir703', 'Mir142hg'], cmap='Reds', ncols=5, vmax='p95')

# %%
sc.pl.umap(adata, color=['Mir124-2hg', 'Mir124a-1hg', '2900055J20Rik', 'Snhg11', 'Syt1', 
                         'Gabra1', 'Cx3cl1', 'Slc12a5', 'Slc17a7', 'Adcy1'], cmap='Reds', ncols=5, vmax='p95')

# %%
from sklearn.covariance import GraphicalLassoCV

# %%
# %%time

model = GraphicalLassoCV(verbose=True, cv=5)
model.fit(adata.X.A)

# %%
import scipy

m = sorted(list(set(adata.var_names[adata.var.highly_variable].tolist()) | set(['Mir9-3hg', 'Mir703', 'Mir124-2hg', 'Mir124a-1hg', 'Mir142hg'])))
ad = adata[:, m]
ad._inplace_subset_obs((ad.X.sum(1)>0).A1)
rho, p = scipy.stats.spearmanr(ad[:, m].X.A)

# %%
m_df = pd.DataFrame(rho, index=m, columns=m)[['Mir9-3hg', 'Mir703', 'Mir124-2hg', 'Mir124a-1hg', 'Mir142hg']]
p_df = pd.DataFrame(p, index=m, columns=m)[['Mir9-3hg', 'Mir703', 'Mir124-2hg', 'Mir124a-1hg', 'Mir142hg']]

# %%
m_df.sort_values('Mir124a-1hg', ascending=False).head(30)

# %%
sc.pl.umap(adata, color=['sample_name', 'predicted_sex', 'leiden'], size=1, ncols=4)

# %% [markdown]
# ## Velocyto

# %%
adata.X = adata.layers['count']

# %%
scv.pp.filter_and_normalize(adata)
adata

# %% [markdown]
# ## Add markers back to matrix

# %%
important_markers = sorted(set(list('''
C1qc Sirt2 Meg3 Pdgfra Gfap Cldn5 C1qb Cnp Gad1 Pcp4 Olig1 Slc1a3 Vwa1 C1qa Mbp Nap1l5 Syt1 Pllp Gstm1 Ctla2a Ctss 
Slc17a7 Htra1 Sparc Hexb Grin2b Atp1a2 Fcer1g Neurod2 Tyrobp Neurod6 Fcrls Satb2 Fezf2
'''.split())))

# %%
hvg_and_markers = sorted(set(list(adata.var_names[adata.var.highly_variable].tolist() + important_markers)))

# %%
adata = adata[:, hvg_and_markers].copy()
adata

# %%
# %%time

scv.pp.moments(adata, n_neighbors=15, n_pcs=50)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
scv.tl.terminal_states(adata)

# %%
sc.pl.umap(adata, color='Cell type', size=2)

# %%
scv.pl.velocity_embedding_grid(adata, basis='umap', color='Cell type', scale=0.3, legend_loc='right margin', alpha=1, size=1)

# %%
#scv.pl.velocity_embedding_stream(adata[~adata.obs['Cell type'].isnull()], basis='umap', color='Cell type', legend_loc = 'right margin')

# %%
sc.pl.umap(adata, color=['root_cells', 'end_points'], cmap='RdYlBu')

# %%
#scv.pl.velocity(adata, var_names=important_markers)
scv.pl.velocity(adata, var_names=['Gad1', 'Neurod2', 'Neurod6'])

# %% [markdown]
# ## Per group velocity

# %%
group_key = 'postnatal_age'
per_group_velocyto = True

# %%
ad.obs.merge(labels[labels.postnatal_age == 21], how='left', on='leiden').set_index(ad.obs.index)

# %%
labels.rename(columns={'sample_leiden': 'leiden'}, inplace=True)

# %%
ad.obs

# %%
for group in (1, 7, 21):
    ad = sc.read(f'adata-{group}-velocyto.h5ad')
    ad.obs.leiden = ad.obs.leiden.astype(str)
    labels.leiden = labels.leiden.astype(str)
    
    del ad.obs['Cell type']
    del ad.obs['postnatal_age_y']
    del ad.obs['postnatal_age_x']
    
    ad.obs = ad.obs.merge(labels[labels.postnatal_age == group], how='left', on=['leiden']).set_index(ad.obs.index)
    display(ad.obs.head())
    ad.write(f'adata-{group}-velocyto.h5ad')

# %%
# %%time

from IPython.core.display import display, HTML
sc.settings.verbosity = 0

individual_ads = {}

if per_group_velocyto:
    ascat = pd.Categorical(adata.obs[group_key])
    for group in tqdm(ascat.categories):
        display(HTML(f'<h2>{group_key}: {group}</h2>'))

        ad = sc.read(f'adata-{group}-velocyto.h5ad')
        ad.obs = ad.obs.merge(labels[labels.postnatal_age == group], how='left', on='leiden').set_index(ad.obs.index)
        sc.pl.umap(ad, color='Cell type')
        
        ad.X = ad.layers['count']
        scv.pp.filter_and_normalize(ad)

        ## Add markers back to matrix
        hvg_and_markers = sorted(set(list(ad.var_names[ad.var.highly_variable].tolist() + important_markers)))
        ad = ad[:, hvg_and_markers].copy()

        scv.pp.moments(ad, n_neighbors=15, n_pcs=50)
        scv.tl.velocity(ad)
        scv.tl.velocity_graph(ad)
        scv.tl.terminal_states(ad)
        
        individual_ads[group] = ad
        
        scv.pl.velocity_embedding_grid(ad, basis='umap', color='Cell type', scale=0.3, legend_loc='right margin', alpha=1, size=1)
        sc.pl.umap(ad, color=['root_cells', 'end_points'], cmap='RdYlBu')
        #scv.pl.velocity(ad, var_names='Fezf2')

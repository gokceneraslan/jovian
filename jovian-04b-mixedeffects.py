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
# # Mixed effect modeling with Scanpy and lme4 R package

# %%
import scanpy as sc
import pandas as pd

# %%
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, r, Formula
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter

def fit_lme(formula, df, nb=False, **fit_kwargs):
    f = Formula(formula)
    
    lmer = importr('lmerTest') # overloads lmer function from lme4 package
    base = importr('base')
    stats = importr('stats')

    with localconverter(ro.default_converter + pandas2ri.converter):
        if not nb:
            fit = lmer.lmer(f, df, **fit_kwargs)
        else:
            fit = r('lme4::glmer.nb')(f, df, **fit_kwargs)
        anova_df = stats.anova(fit)

    # no automatic pandas conversion here, keep as FloatMatrix
    coef_df = base.summary(fit).rx2('coefficients')

    with localconverter(ro.default_converter + pandas2ri.converter):
        coef_df = r['as.data.frame'](coef_df)

    return coef_df, anova_df


# %%
from tqdm.auto import tqdm

def fit_lme_adata(adata, genes, formula, obs_features, use_raw=False):
    coef_df = {}
    anova_df = {}
    covariates = adata.obs[obs_features]
    
    for gene in tqdm(genes):
        gene_vec = adata[:, gene].X if not use_raw else adata.raw[:, gene].X
        df = pd.concat([covariates,
                        pd.DataFrame(gene_vec,
                                     index=adata.obs.index, 
                                     columns=['gene'])], axis=1)

        coefs, anova = fit_lme(formula, df)
        coef_df[gene] = coefs
        anova_df[gene] = anova
        
    coef_df = pd.concat([df.assign(gene=gene) for gene, df in coef_df.items()], axis=0)
    coef_df = coef_df.reset_index().rename(columns={'index': 'fixed_effect'})

    anova_df = pd.concat([df.assign(gene=gene) for gene, df in anova_df.items()], axis=0)
    anova_df = anova_df.reset_index().rename(columns={'index': 'fixed_effect'})
   
    return coef_df, anova_df


# %%
adata = sc.datasets.paul15()
adata.obs['somecoef'] = adata[:, 'Cst3'].X
adata.obs.head()

# %%
c, a = fit_lme_adata(adata, ['Gata1', 'Gata2'], 'gene ~ somecoef + (1|paul15_clusters)', ['paul15_clusters', 'somecoef'])

# %%
c

# %%
a

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import scvelo as scv

scv.logging.print_version()
sc.settings.autoshow = False
sc.set_figure_params(dpi_save=300, dpi=80)

def scanpy_bbknn_pipeline(adata, neighbors=3, PCs=30, HVGs=2000, regress=True, regress_list=['total_counts','pct_counts_mt'], bbknn=True, bbknn_batch_key='batch', n_neighbors=10):
    adata.raw = adata
    
    sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.highly_variable_genes(adata, n_top_genes=HVGs) 
    adata = adata[:, adata.var.highly_variable]
    if regress:
        sc.pp.regress_out(adata, regress_list, n_jobs=n_jobs)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=True, n_comps=PCs)
    sc.external.pp.bbknn(adata, n_pcs=PCs, neighbors_within_batch=neighbors, batch_key=bbknn_batch_key)
    sc.tl.umap(adata)
    
    return(adata)

adata_tmp = sc.read_h5ad('tmp.h5ad')
adata = adata_tmp[adata_tmp.obs['subCluster'].isin(['TE12-CD8-Tem-CXCR4', 'TE13-CD8-Tem-IFNG', 'TE14-CD8-Tex-CXCL13'])]
adata = adata[adata.obs['Tissue']=='Rectum_T']
adata = scanpy_bbknn_pipeline(adata)
scv.pp.filter_and_normalize(adata,
                            min_shared_cells=None,
                            n_top_genes = 2000
                           )
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(adata, n_jobs=40)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata, n_jobs=40)
scv.pl.velocity_embedding_stream(adata, basis='umap', color='subCluster', arrow_color="black", legend_loc='right margin', legend_fontsize=5, density=1.0, save='RNA_velocity_stream.pdf')
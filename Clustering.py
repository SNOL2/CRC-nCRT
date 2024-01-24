import sys, os, re, gc, pickle, math, time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from collections import Counter


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


adata = sc.read_h5ad('tmp.h5ad')
adata = scanpy_bbknn_pipeline(adata)
sc.tl.leiden(adata)
sc.tl.rank_genes_groups(adata, 'leiden')
sc.pl.umap(adata, color=[ "EPCAM","PECAM1","DCN","PTPRC","CPA3","CD14","LYZ","CD3D","CD3G","CD4","CD8A","NKG7","NCAM1","MS4A1","MZB1",'leiden', 'Ident', 'Tissue'],legend_loc="on data" )
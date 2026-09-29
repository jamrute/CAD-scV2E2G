#!/usr/bin/env python
"""
03_perturbseq/06_cnmf/02_make_h5ad.py
Purpose : Build AnnData from exported counts; drop all-zero genes and genes
          expressed in < 5 cells (Methods).
Status  : REFERENCE IMPLEMENTATION
"""
import argparse
import numpy as np
import pandas as pd
import scipy.io
import anndata as ad

ap = argparse.ArgumentParser()
ap.add_argument("--indir", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--min-cells", type=int, default=5)
a = ap.parse_args()

X = scipy.io.mmread(f"{a.indir}/matrix.mtx").T.tocsr()          # cells x genes
genes = pd.read_csv(f"{a.indir}/genes.tsv", header=None)[0].values
cells = pd.read_csv(f"{a.indir}/cells.tsv", header=None)[0].values
obs = pd.read_csv(f"{a.indir}/obs.tsv", sep="\t").set_index("cell").loc[cells]

adata = ad.AnnData(X=X.astype(np.float32), obs=obs, var=pd.DataFrame(index=genes))
n_cells = np.asarray((adata.X > 0).sum(axis=0)).ravel()
adata = adata[:, n_cells >= a.min_cells].copy()
print(f"{adata.n_obs} cells x {adata.n_vars} genes retained")
adata.write_h5ad(a.out)

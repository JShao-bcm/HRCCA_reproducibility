import sys
import scanpy as sc
import anndata
import pandas as pd
import os
import numpy as np
from pathlib import Path

#mamba activate scglue_jj/

main="/dfs3b/ruic20_lab/jianms1/Single_cell/scRNA_snRNA_1-7-2025/snRNA_NoFetal_clean/annotation/"
#dir1="/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/genexp_sample_raw_atlas"
#dir1=f"{main}/7_DEG/genexp_sample_raw_atlas"
dir1=f"{main}/Aging_LMM/Jun_pipeline/1_genexp_donor_raw_atlas"

Path(dir1).mkdir(parents=True,exist_ok=True)

#path=f"{main}/5_refine_major/scvi/major/ON_ONH_clean_new.h5ad"
#path=f"{main}/5_refine_major/scvi/major/major_ON_ONH_new_final.h5ad"
path=f"{main}/snRNA_clean_all_rawcounts_final_mt_filtered.h5ad"
adata=sc.read_h5ad(path)

snRNA_majorclass_new = {"Macrophages": "Macrophages",
                        "Plasma cells": "Plasma_cells",
                        "B cells": "B_cells",
                        "Monocytes": "Monocytes",
                        "T cells": "T_cells",
                        "NK cells": "NK_cells",
                        "Melanocytes": "Melanocytes",
                        "RPE cells": "RPE_cells",
                        "Schwann cells": "Schwann_cells",
                        "Fibroblasts": "Fibroblasts",
                        "Mural cells": "Mural_cells",
                        "Endothelial cells": "Endothelial_cells",
                        "Mast cells": "Mast_cells",
                        "DC": "DC",
                        "pDC": "pDC",
                        "ILC3": "ILC3",
                        "Megakaryocytes": "Megakaryocytes"
                       }

adata.obs["snRNA_majorclass_new"] = adata.obs.snRNA_majorclass.map(snRNA_majorclass_new)
#adata.obs["donor"]=adata.obs["donor"].replace("Chen_17D013","Chen_D017_13")

def rename_gene(adata):
	for i in range(0, adata.var.features.shape[0]):
		if adata.var.features[i] in gene_dict:
			adata.var.features[i] = gene_dict[adata.var.features[i]]
	return adata

def grouped_obs_cpm(adata, group_key, layer=None, gene_symbols=None):
	if layer is not None:
		getX= lambda x: x.layers[layer]
	else:
		getX = lambda x: x.raw.X
	if gene_symbols is not None:
		new_idx = adata.var[idx]
	else:
		new_idx = adata.var_names
	grouped = adata.obs.groupby(group_key)
	out = pd.DataFrame(
		np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
		columns=list(grouped.groups.keys()),
		index=adata.var_names
	)
	for group, idx in grouped.indices.items():
		X = getX(adata[idx])
		total_count=X.sum(dtype=np.float64)
#		factor=1000000.0
#		out[group] = np.ravel(X.sum(axis=0, dtype=np.float64))/total_count*factor
		out[group] = np.ravel(X.sum(axis=0, dtype=np.float64))
	return(out)

#celltype=list(adata.obs["snRNA_majorclass_new"].unique())
celltype=['Macrophages','Plasma_cells', 'B_cells', 'Monocytes', 'T_cells', 'NK_cells', 'Melanocytes', 'RPE_cells', 'Schwann_cells', 'Fibroblasts', 'Mural_cells', 'Endothelial_cells']

adata=adata[adata.obs.snRNA_majorclass_new.isin(celltype)].copy()

cellnum={}

list1=f"{main}/Aging_LMM/Jun_pipeline/snRNA_donor_majorclass_num_atlas_cell"

out=open(list1,"w")

for donor in adata.obs["Donor"].value_counts().keys() :
		adata_sub = adata[adata.obs["Donor"]==donor]
		adata_sub.layers["raw"] = adata_sub.X
		genexp=grouped_obs_cpm(adata=adata_sub, group_key="snRNA_majorclass_new", layer="raw")
		fileout=f'{dir1}/{donor}.txt.gz'
		genexp.to_csv(fileout,sep="\t")
		for cell in celltype :
			if cell in adata_sub.obs["snRNA_majorclass_new"].value_counts():
				if donor not in cellnum:
					cellnum[donor]={}
				cellnum[donor][cell] = adata_sub.obs["snRNA_majorclass_new"].value_counts()[cell]
				num_tmp=adata_sub.obs["snRNA_majorclass_new"].value_counts()[cell]
				out.write(f'{donor}\t{cell}\t{num_tmp}\n')
			else:
				if donor not in cellnum:
					cellnum[donor]={}
				cellnum[donor][cell] = 0
				num_tmp=0
				out.write(f'{donor}\t{cell}\t{num_tmp}\n')

out.close()

df=pd.DataFrame.from_dict(cellnum)

list=f"{main}/Aging_LMM/Jun_pipeline/snRNA_donor_majorclass_num_atlas"

df.to_csv(f'{list}')


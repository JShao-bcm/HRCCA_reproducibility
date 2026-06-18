import sys
import os
sys.path.insert(0, os.path.abspath("/pub/jianms1/software/NSForest/"))
sys.path.insert(0, os.path.abspath("/pub/jianms1/software/NSForest/nsforest/nsforesting"))
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import nsforest as ns
from nsforest import utils
from nsforest import preprocessing as pp
from nsforest import nsforesting
from nsforest import evaluating as ev
from nsforest import plotting as pl

in_dir=sys.argv[1] #/dfs3b/ruic20_lab/jianms1/Single_cell/Stromal_cells/
bn=sys.argv[2] #base name of input file: eg. Stromal_cells_norm
cluster_header = sys.argv[3]#'subclass_V1' 

file=f'{in_dir}/annotation/{bn}.h5ad'

adata = sc.read_h5ad(file)
adata

output_folder = f"{in_dir}/annotation/NSForest/{cluster_header}/"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
   
os.chdir(output_folder)

adata.obs_names
adata.var_names

adata.obs[cluster_header]=adata.obs[cluster_header].astype('category')
adata.obs[cluster_header].value_counts()

ns.pp.dendrogram(adata, cluster_header, save_plot = False, output_folder = output_folder, outputfilename_suffix = cluster_header)
adata = ns.pp.prep_medians(adata, cluster_header)
adata

adata = ns.pp.prep_binary_scores(adata, cluster_header)
adata

plt.clf()
filename = output_folder + cluster_header + '_medians.png'
print(f"Saving median distributions as...\n{filename}")
a = plt.figure(figsize = (6, 4))
a = plt.hist(adata.varm["medians_" + cluster_header].unstack(), bins = 100)
a = plt.title(f'{file.split("/")[-1].replace(".h5ad", "")}: {"medians_" + cluster_header} histogram')
a = plt.xlabel("medians_" + cluster_header)
a = plt.yscale("log")
a = plt.savefig(filename, bbox_inches='tight')
#plt.show()

plt.clf()
filename = output_folder + cluster_header + '_binary_scores.png'
print(f"Saving binary_score distributions as...\n{filename}")
a = plt.figure(figsize = (6, 4))
a = plt.hist(adata.varm["binary_scores_" + cluster_header].unstack(), bins = 100)
a = plt.title(f'{file.split("/")[-1].replace(".h5ad", "")}: {"binary_scores_" + cluster_header} histogram')
a = plt.xlabel("binary_scores_" + cluster_header)
a = plt.yscale("log")
a = plt.savefig(filename, bbox_inches='tight')
#plt.show()

filename = file.replace(".h5ad", "_preprocessed_{res}.h5ad")
print(f"Saving new anndata object as...\n{filename}")
adata.write_h5ad(filename)

outputfilename_prefix = cluster_header
results = nsforesting.NSForest(adata, cluster_header, output_folder = output_folder, gene_selection = 'BinaryFirst_moderate',outputfilename_prefix = outputfilename_prefix)

#results.to_csv(f'immune_majorclass.csv',sep=',',index=False)

ns.pl.boxplot(results, "f_score",save=True,output_folder='./',outputfilename_prefix='f_score.boxplot')
ns.pl.boxplot(results, "PPV",save=True,output_folder='./',outputfilename_prefix='ppv.boxplot')
ns.pl.boxplot(results, "recall",save=True,output_folder='./',outputfilename_prefix='recall.boxplot')
ns.pl.boxplot(results, "onTarget",save=True,output_folder='./',outputfilename_prefix='onTarget.boxplot')
ns.pl.scatter_w_clusterSize(results, "f_score",save=True,output_folder='./',outputfilename_prefix='f_score.scatter')
ns.pl.scatter_w_clusterSize(results, "PPV",save=True,output_folder='./',outputfilename_prefix='ppv.scatter')
ns.pl.scatter_w_clusterSize(results, "recall",save=True,output_folder='./',outputfilename_prefix='recall.scatter')
ns.pl.scatter_w_clusterSize(results, "onTarget",save=True,output_folder='./',outputfilename_prefix='onTarget.scatter')

to_plot = results.copy()
dendrogram = [] # custom dendrogram order
dendrogram = list(adata.uns["dendrogram_" + cluster_header]["categories_ordered"])
to_plot["clusterName"] = to_plot["clusterName"].astype("category")
to_plot["clusterName"] = to_plot["clusterName"].cat.set_categories(dendrogram)
to_plot = to_plot.sort_values("clusterName")
to_plot = to_plot.rename(columns = {"NSForest_markers": "markers"})
to_plot.head()

markers_dict = dict(zip(to_plot["clusterName"], to_plot["markers"]))
markers_dict

ns.pl.dotplot(adata, markers_dict, cluster_header, dendrogram = dendrogram, save = True, output_folder = output_folder, outputfilename_suffix = outputfilename_prefix)
ns.pl.stackedviolin(adata, markers_dict, cluster_header, dendrogram = dendrogram, save = True, output_folder = output_folder, outputfilename_suffix = outputfilename_prefix)
ns.pl.matrixplot(adata, markers_dict, cluster_header, dendrogram = dendrogram, save = True, output_folder = output_folder, outputfilename_suffix = outputfilename_prefix)


#args=commandArgs(trailingOnly=T)
#ct=read.table(args[1],header=F)$V1
ct=c('Macrophages','Plasma_cells', 'B_cells', 'Monocytes', 'T_cells', 'NK_cells', 'Melanocytes', 'RPE_cells', 'Schwann_cells', 'Fibroblasts', 'Mural_cells', 'Endothelial_cells')
#dir="/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/genexp_sample_raw_celltype_atlas/"
dir="2_genexp_donor_count_matrix" #args[2]

#cell=c("Oligodendrocyte","Astrocyte","Fibroblast","Rod","Microglia","Endothelial_cell","BC","Oligodendrocyte_precursor_cell","Mural_cell","MG","AC","RPE","Macrophage","Melanocyte","HC","Cone","RGC","Schwann_cell")
#dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas"
for(c in ct){
file=paste0(dir,"/exp_",c,"_Norm_all" )

data=read.table(file,header=T)
cor=cor(data)

out=paste0(dir,"/exp_",c,"_cor","_Norm_all" )


write.table(cor,file=out,quote=F,sep="\t")
}

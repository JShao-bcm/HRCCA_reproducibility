library('edgeR')

setwd("/dfs3b/ruic20_lab/jianms1/Single_cell/scRNA_snRNA_1-7-2025/snRNA_NoFetal_clean/annotation/Aging_LMM/Jun_pipeline")

#args=commandArgs(trailingOnly=T)
#ct=read.table(args[1],header=F)$V1
ct=c('Macrophages','Plasma_cells', 'B_cells', 'Monocytes', 'T_cells', 'NK_cells', 'Melanocytes', 'RPE_cells', 'Schwann_cells', 'Fibroblasts', 'Mural_cells', 'Endothelial_cells')
seq="all"
#dir="/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/genexp_sample_raw_celltype_atlas/"
dir="2_genexp_donor_count_matrix" #args[2] #"/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas"

for(cell in ct){
exp=read.table(paste0(dir,"/exp_",cell),header=T)
fListNames=colnames(exp)

fListNames=gsub("X","",fListNames)
fListNames=gsub("\\.","\\-",fListNames)

exp1=exp
countMatrix=data.matrix(exp1)
t=quantile(rowMeans(cpm(countMatrix)),seq(0,1,0.05))
######isexpr=rowMeans(cpm(countMatrix)) >= t[15]

######geneExpr = DGEList( countMatrix[isexpr,] )
geneExpr = DGEList( countMatrix )

geneExpr = calcNormFactors( geneExpr )

CPM <- cpm(geneExpr, prior.count=0, log=F)


write.table(CPM,file=paste0(dir,"/exp_",cell,"_Norm_",seq),sep="\t",quote=F)


logCPM <- cpm(geneExpr, prior.count=1, log=TRUE)

write.table(logCPM,file=paste0(dir,"/exp_",cell,"_logNorm_",seq),sep="\t",quote=F)

}


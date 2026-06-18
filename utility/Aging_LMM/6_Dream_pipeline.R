library('variancePartition')
library('edgeR')
library('BiocParallel')
#library("DESeq2")
library("ggplot2")
#library( "gplots" )
library( "RColorBrewer" )
#library( "genefilter" )
#library("ggrepel")
library(tidyr)
#library(peer)
library(qvalue)

flt_count = function(ct_exp, cpm_exp,min_ct, min_cpm,min_samples)
{
exp_ct=as.matrix(ct_exp)
exp_ct_raw=as.matrix(ct_exp)
exp_cpm=as.matrix(cpm_exp)
exp_ct[exp_ct < min_ct]=0
exp_ct[exp_ct >= min_ct]=1
exp_cpm[exp_cpm < min_cpm]=0
exp_cpm[exp_cpm >= min_cpm]=1
exp_f = exp_ct_raw[(rowSums(exp_ct)>= (min_samples * ncol(exp_ct)))&(rowSums(exp_cpm)>= (min_samples * ncol(exp_cpm))), ]
return(exp_f)
}



############
args <- commandArgs(trailingOnly = TRUE)

cell=args[1] #cell type name, eg Endothelial_cells
dirInName=args[2] #maybe output of cell type dreamea results
fn=args[3] #metadata, /dfs3b/ruic20_lab/jianms1/Single_cell/scRNA_snRNA_1-7-2025/snRNA_NoFetal_clean/annotation/Aging_LMM/scRNA_snRNA_Demographic_Ethnicity.csv
rml=read.table(args[4],header=T) #filter sample file, /dfs3b/ruic20_lab/jianms1/Single_cell/scRNA_snRNA_1-7-2025/snRNA_NoFetal_clean/annotation/Aging_LMM/Jun_pipeline/2_genexp_donor_count_matrix/exp_{celltype}_cor_Norm_all_flt
dirIn=args[5] #output dir, eg. /dfs3b/ruic20_lab/jianms1/Single_cell/scRNA_snRNA_1-7-2025/snRNA_NoFetal_clean/annotation/Aging_LMM/Jun_pipeline/3_Dream_majorclass/Dream_donor_DEG_{celltype}
indir=args[6] #dir to gene exp rawcount, /dfs3b/ruic20_lab/jianms1/Single_cell/scRNA_snRNA_1-7-2025/snRNA_NoFetal_clean/annotation/Aging_LMM/Jun_pipeline/2_genexp_donor_count_matrix

dir.create(dirIn,recursive = TRUE)
setwd(dirIn)
dirIn <- getwd()
exp=read.table(paste0(indir,"/exp_",cell),header=T)
file_info=read.table(paste0(fn),header=T,sep=",")
#file_info=file_info[(file_info$age_year!="Unk")&(file_info$gender!="Unk")&(!(file_info$sampleid %in% rml$V1))&(!(file_info$sampleid %in% mismatch_sample)),]
file_info=file_info[(file_info$Age!="Unk")&(file_info$Gender!="Unk")&(!(file_info$Donor %in% rml$V1)),]#&(!(file_info$Donor %in% mismatch_sample)),]

fListNames=colnames(exp)
fListNames=gsub("X","",fListNames)
fListNames=gsub("\\.","\\-",fListNames)

m=match(file_info$Donor,fListNames,nomatch=0)

exp=exp[,m]
fListNames=colnames(exp)
fListNames=gsub("X","",fListNames)
fListNames=gsub("\\.","\\-",fListNames)
colnames(exp)=fListNames
m=match(fListNames,file_info$Donor,nomatch=0)

metadata=file_info[m,]
metadata$age=as.numeric(metadata$Age)

metadata$race=factor(metadata$Inferred,levels=c("Black","White","Asian","Hispanic"))
metadata$gender=factor(metadata$Gender,levels=c("Male","Female"))

#metadata$tissue=factor(metadata$tissue,levels=c("ON","ONH"))

rownames(metadata)=metadata$Donor


#isexpr=rowMeans(cpm(exp)) > 0.5 #t[15]  #rowMeans >=1
#isexpr=rowMeans(cpm(exp)) >= 1 #t[15]  #rowMeans >=1
isexpr=rowMeans(cpm(exp)) >= 5 #t[15]  #rowMeans >=1

countMatrix=exp[isexpr,]
#countMatrix=exp
#geneExpr = DGEList( countMatrix[isexpr,] )

geneExpr = DGEList( countMatrix )
geneExpr = calcNormFactors( geneExpr )

#########

metadata$libsize=floor(log(geneExpr$samples[metadata$Donor,]$lib.size,base=10))

library(dplyr)
library(stringr)

metadata <- metadata %>%
  mutate(source = case_when(
    str_detect(Donor, "BCM") ~ "BCM",
    str_detect(Donor, "MMD") ~ "MMD",
    TRUE                     ~ "Sanes"  # Everything else is categorized as Sanes
  ))

#########

param = SnowParam(4, "SOCK", progressbar=TRUE)

form <- ~ age  + gender +  race +  (1|source) ### all_new_all_mac_clean

vobjDream = voomWithDreamWeights( geneExpr, form, metadata, BPPARAM=param ) 
fitmm = dream( vobjDream, form, metadata)

#fitmm = dream( vobjDream, form, metadata, L)
fitmm = eBayes(fitmm)

head(fitmm$design, 3)

saveRDS(vobjDream,file=paste0(cell,"_vobjDream_cpm5_no_libsize.rds"))
saveRDS(fitmm,file=paste0(cell,"_fitmm_cpm5_no_libsize.rds"))

res=topTable( fitmm, coef='age', number=Inf )
res$qval=qvalue(res$P.Value)$qvalues
#write.table(res,file=paste0(cell,"_DEG_age"),sep="\t",quote=F)
write.table(res,file=paste0(cell,"_DEG_age_cpm5_no_libsize"),sep="\t",quote=F)

res=topTable( fitmm, coef='genderFemale', number=Inf )
res$qval=qvalue(res$P.Value)$qvalues
#write.table(res,file=paste0(cell,"_DEG_gender"),sep="\t",quote=F)
write.table(res,file=paste0(cell,"_DEG_gender_cpm5_no_libsize"),sep="\t",quote=F)

res=topTable( fitmm, coef='raceWhite', number=Inf )
res$qval=qvalue(res$P.Value)$qvalues
##write.table(res,file=paste0(cell,"_DEG_raceBlack"),sep="\t",quote=F)
write.table(res,file=paste0(cell,"_DEG_raceWhite_cpm5_no_libsize"),sep="\t",quote=F)

res=topTable( fitmm, coef='raceHispanic', number=Inf )
res$qval=qvalue(res$P.Value)$qvalues
#write.table(res,file=paste0(cell,"_DEG_raceHispanic"),sep="\t",quote=F)

write.table(res,file=paste0(cell,"_DEG_raceHispanic_cpm5_no_libsize"),sep="\t",quote=F)

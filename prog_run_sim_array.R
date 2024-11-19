library(tidyverse)
library(reshape2)
library(lme4)
library(SummarizedExperiment)
library(DESeq2)
library(edgeR)
library(parallel)
library(doParallel)

load("sims_N100_m2_s1_rtruncnorm13.rda")
ind_run=as.integer(commandArgs(trailingOnly = TRUE)[1]) #can be run in parallel using array on BH
ncores=as.integer(commandArgs(trailingOnly = TRUE)[2])
print(ind_run)
print(ncores)

#start cluster to run miRNA in parallel
#ncores=8 #match to what was requested on BH
if (ncores>1){
  cl=makeCluster(ncores)
  registerDoParallel(cl)
}
#startTime=Sys.time()

#define groups to compare
col_group_in = c(rep("A", 19), rep("B",20))

source('miRglmm.R')
source('miRglm.R')

#fit miRglmm full and reduced models
fits = miRglmm(sims[[ind_run]]$sim_se, col_group=col_group_in, ncores = ncores)

#fit miRglmm full and reduced models
fits[["miRglmm poisson"]] = miRglmm(sims[[ind_run]]$sim_se, col_group=col_group_in, ncores = ncores, family="poisson")

#aggregate data to miRNAs and fit miRglmnb
miRNA_counts = t(apply(assay(sims[[ind_run]]$sim_se), 2, function(x) by(x, rowData(sims[[ind_run]]$sim_se)$miRNA, sum)))
fits[["miRglmnb"]]= miRglm(miRNA_counts, col_group=col_group_in, ncores = ncores)
fits[["miRglmpois"]]= miRglm(miRNA_counts, col_group=col_group_in, ncores = ncores, family="poisson")


if (ncores>1){
  stopCluster(cl)
}
# run wilcoxon 

fits[["wilcoxon"]]=data.frame("wilcoxp"=apply(t(miRNA_counts), 1, function(x) wilcox.test(x[which(col_group_in=="A")], x[which(col_group_in=="B")])$p.value))


#run DESeq2 
panelB_raw=SummarizedExperiment(assays=list(t(miRNA_counts)), rowData=colnames(miRNA_counts), colData=col_group_in)
names(colData(panelB_raw))=c("col_group")
ddsSE=DESeqDataSet(panelB_raw, design=~col_group)
fits[["DESeq2"]]=DESeq(ddsSE)

#run edgeR
edgeR_set=DGEList(counts=t(miRNA_counts), group=col_group_in)
design=model.matrix(~col_group_in)
edgeR_set=estimateDisp(edgeR_set, design)
et=glmQLFit(edgeR_set, design)
fits[["edgeR"]]=glmQLFTest(et)


#run limma-voom
limvoom_set=DGEList(counts=t(miRNA_counts))
y=voom(limvoom_set, design) #same design as edgeR
limvoom_fit=lmFit(y,design)
fits[["limvoom"]]=eBayes(limvoom_fit, trend=TRUE)


#will save file containing one list of 6 model fits for each simulated dataset
save(fits, file=paste0("sim results/sims_N100_m2_s1_rtruncnorm13_results", ind_run, ".rda"))


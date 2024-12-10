library(tidyverse)
library(reshape2)
library(lme4)
library(SummarizedExperiment)
library(DESeq2)
library(edgeR)

#load truth data
ratio_pool=read.csv(file="ERCC files/FINAL_Ratiometric_SynthA_and_SynthB-1_fixlast3.csv", sep="\t")
ratio_pool$A=as.numeric(gsub("x","",ratio_pool$X.SAMPLE.A.))
ratio_pool$B=as.numeric(gsub("x","",ratio_pool$X.SAMPLE.B.))
ratio_pool$ratio=ratio_pool$B/ratio_pool$A
ratio_pool$true_logFC=log(ratio_pool$ratio)
out=ratio_pool["true_logFC"]
rownames(out)=ratio_pool$ratio.seqID



#load processed results
load('ERCC/allfilters_processed_results.rda')


beta_hat=results[["beta_hat"]]
beta_hat=transform(merge(beta_hat, out, by='row.names'), row.names=Row.names, Row.names=NULL)
beta_hat=beta_hat[,c("filter..1","miRglmm_poiss", "miRglmnb","miRglmpois", "DESeq2","edgeR","limmavoom", "true_logFC")]
colnames(beta_hat)=c("miRglmm-NB","miRglmm-Poisson", "NB GLM","Poisson GLM", "DESeq2","edgeR","limma-voom", "True LogFC")


diff_squared=(beta_hat[,!(colnames(beta_hat)=="True LogFC")]-beta_hat[,"True LogFC"])^2
MSE=data.frame(MSE=colMeans(diff_squared, na.rm=TRUE))
MSE$MSE=MSE$MSE*1000


CI_LL=results[["CI_LL"]]
CI_LL=transform(merge(CI_LL, out, by='row.names'), row.names=Row.names, Row.names=NULL)
CI_LL=CI_LL[,c("filter..1","miRglmm_poiss", "miRglmnb","miRglmpois", "DESeq2","limmavoom", "true_logFC")]
colnames(CI_LL)=c("miRglmm-NB","miRglmm-Poisson", "NB GLM","Poisson GLM", "DESeq2","limma-voom", "True LogFC")
gr_LL=CI_LL[,!(colnames(CI_LL)=="True LogFC")]<=CI_LL[,"True LogFC"]


CI_UL=results[["CI_UL"]]
CI_UL=transform(merge(CI_UL, out, by='row.names'), row.names=Row.names, Row.names=NULL)
CI_UL=CI_UL[,c("filter..1","miRglmm_poiss", "miRglmnb","miRglmpois", "DESeq2","limmavoom", "true_logFC")]
colnames(CI_UL)=c("miRglmm-NB","miRglmm-Poisson", "NB GLM","Poisson GLM", "DESeq2","limma-voom", "True LogFC")
lt_UL=CI_UL[,!(colnames(CI_UL)=="True LogFC")]>=CI_UL[,"True LogFC"]

cov_prob=data.frame("coverage probability"=colMeans(gr_LL& lt_UL, na.rm=TRUE))

nullvar_mat=data.frame("null variance"=apply(beta_hat[beta_hat$`True LogFC`==0,],2,var, na.rm=TRUE))
nullvar_mat=nullvar_mat*1000
rownames(nullvar_mat)=c("miRglmm-NB","miRglmm-Poisson", "NB GLM","Poisson GLM", "DESeq2","edgeR", "limma-voom", "True LogFC")



DEvar_mat=data.frame("null variance"=apply(rbind(-1*beta_hat[beta_hat$`True LogFC`<0,], beta_hat[beta_hat$`True LogFC`>0,]),2,var, na.rm=TRUE))
DEvar_mat=DEvar_mat*1000
rownames(DEvar_mat)=c("miRglmm-NB","miRglmm-Poisson", "NB GLM","Poisson GLM", "DESeq2","edgeR", "limma-voom", "True LogFC")

pvals=results[["pvals"]]
all(rownames(pvals)==rownames(beta_hat))
pvals=pvals[,c("filter..1","miRglmm_poiss", "miRglmnb","miRglmpois", "DESeq2","edgeR","limmavoom")]
colnames(pvals)=c("miRglmm-NB","miRglmm-Poisson", "NB GLM","Poisson GLM", "DESeq2","edgeR","limma-voom")


FDR=data.frame(apply(pvals, 2, function(x) p.adjust(x, method="fdr")))
colnames(FDR)=colnames(pvals)
all(rownames(FDR)==rownames(beta_hat))

is_sig=as.data.frame(FDR<0.05)
TPR=data.frame("TPR"=colSums(is_sig[which(beta_hat$`True LogFC`!=0),], na.rm=TRUE)/colSums(!is.na(is_sig[which(beta_hat$`True LogFC`!=0),])))
rownames(TPR)=c("miRglmm-NB","miRglmm-Poisson", "NB GLM","Poisson GLM", "DESeq2","edgeR", "limma-voom")

not_sig=as.data.frame(FDR>0.05)
TNR=data.frame("TNR"=colSums(not_sig[which(beta_hat$`True LogFC`==0),], na.rm=TRUE)/colSums(!is.na(not_sig[which(beta_hat$`True LogFC`==0),])))
rownames(TNR)=c("miRglmm-NB","miRglmm-Poisson", "NB GLM","Poisson GLM", "DESeq2","edgeR", "limma-voom")

FC_vec=sort(unique(log(exp(abs(beta_hat$`True LogFC`)))))
FC_vec=FC_vec[FC_vec>0]
FC_vec=FC_vec[1:7]
for (ind in seq(1, length(FC_vec))){
  FC_in=exp(FC_vec[ind])
  idx=which(exp(abs(beta_hat$`True LogFC`))==FC_in)
  beta_hat_sub=beta_hat[idx,]
  if (ind==1){
  TPR_byFC=data.frame("TPR"=colSums(is_sig[idx,], na.rm=TRUE)/colSums(!is.na(is_sig[idx,])))
  var_byFC=data.frame(apply(rbind(-1*beta_hat_sub[beta_hat_sub$`True LogFC`<0,], beta_hat_sub[beta_hat_sub$`True LogFC`>0,]),2,var)*1000)
  } else {
    TPR_byFC=cbind(TPR_byFC, data.frame("TPR"=colSums(is_sig[idx,], na.rm=TRUE)/colSums(!is.na(is_sig[idx,]))))
    var_byFC=cbind(var_byFC, data.frame(apply(rbind(-1*beta_hat_sub[beta_hat_sub$`True LogFC`<0,], beta_hat_sub[beta_hat_sub$`True LogFC`>0,]),2,var)*1000))
  }
}
colnames(TPR_byFC)=exp(FC_vec)
colnames(var_byFC)=exp(FC_vec)
var_byFC=var_byFC[1:7,]


table_out=cbind(MSE,  "coverage probability"=cov_prob[match(rownames(MSE), rownames(cov_prob)),1])
table_out=cbind(table_out,  "null variance"=nullvar_mat[match(rownames(table_out), rownames(nullvar_mat)),1])
table_out=cbind(table_out,  "DE variance"=DEvar_mat[match(rownames(table_out), rownames(DEvar_mat)),1])
table_out=cbind(table_out,  "TPR"=TPR[match(rownames(table_out), rownames(TPR)),1])
table_out=cbind(table_out,  "TNR"=TNR[match(rownames(table_out), rownames(TNR)),1])


beta_hat$isDE=0
idx=which(beta_hat$`True LogFC`!=0)
beta_hat$isDE[idx]=1
AUC_vec=0
for (ind in seq(1,7)){
  roc_obj=roc(beta_hat$isDE,FDR[,ind])
  AUC_vec[ind]=auc(roc_obj)
}
AUC_vec=data.frame(AUC_vec)
rownames(AUC_vec)=colnames(FDR)
table_out=cbind(table_out,  "AUC"=AUC_vec[match(rownames(table_out), rownames(AUC_vec)),1])

write.csv(table_out, file="figures/Table2.csv")
write.csv(TPR_byFC, file="figures/SupTable27_TPRbyFC.csv")
write.csv(var_byFC, file="figures/SupTable28_varbyFC.csv")
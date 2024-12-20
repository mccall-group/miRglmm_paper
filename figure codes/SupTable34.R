library(ggplot2)

load(file='study89_exact/results/filter_neg1_processed_results.rda')

beta_df=data.frame(results[["beta_hat"]])
beta_df=beta_df[, -c(2,4)]
beta_df[, 1:5]=round(beta_df[, 1:5], 2)
CI_LL_df=data.frame(results[["CI_LL"]])[, -c(2,4)]
CI_LL_df[, 1:4]=round(CI_LL_df[, 1:4],2)
CI_UL_df=data.frame(results[["CI_UL"]])[, -c(2,4)]
CI_UL_df[, 1:4]=round(CI_UL_df[, 1:4],2)
for (ind in seq(1, length(beta_df[1,])-1)){
  method_in=colnames(beta_df)[ind]
  beta_in=beta_df[, method_in]
  if(method_in != "edgeR"){
  LL_in=CI_LL_df[,method_in]
  UL_in=CI_UL_df[,method_in]
  }
  if (ind==1){
  est_comb=paste(beta_in, paste0("(", LL_in, ","), paste0(UL_in, ")"))
  } else if (ind==4){
    est_comb=cbind(est_comb, paste(beta_in, paste0("(", NA, ","), paste0(NA, ")")))
  }
  else{
    est_comb=cbind(est_comb, paste(beta_in, paste0("(", LL_in, ","), paste0(UL_in, ")")))
  }
}
est_comb=data.frame(est_comb)
est_comb$contrast=beta_df$contrast
colnames(est_comb)=colnames(beta_df)
rownames(est_comb)=rownames(beta_df)



p_df=data.frame(results[["pvals"]])
p_df=p_df[, -c(2,4)]
padj_df=matrix(, nrow=dim(p_df)[1], ncol=dim(p_df)[2]-1)
for (ind in seq(1,dim(p_df)[2]-1)){
  padj_df[which(p_df$contrast=="CD4 T lymphocyte vs B lymphocyte"),ind]=p.adjust(p_df[which(p_df$contrast=="CD4 T lymphocyte vs B lymphocyte"),ind], method="BH")
  padj_df[which(p_df$contrast=="CD8 T lymphocyte vs B lymphocyte"),ind]=p.adjust(p_df[which(p_df$contrast=="CD8 T lymphocyte vs B lymphocyte"),ind], method="BH")
  padj_df[which(p_df$contrast=="Monocyte vs B lymphocyte"),ind]=p.adjust(p_df[which(p_df$contrast=="Monocyte vs B lymphocyte"),ind], method="BH")
  padj_df[which(p_df$contrast=="Natural Killer vs B lymphocyte"),ind]=p.adjust(p_df[which(p_df$contrast=="Natural Killer vs B lymphocyte"),ind], method="BH")
}
sig_df=data.frame(padj_df<0.05)
sig_df$contrast=p_df$contrast
colnames(sig_df)=colnames(p_df)
rownames(sig_df)=rownames(data.frame(results[["pvals"]]))


all(rownames(est_comb)==rownames(sig_df))
idx_keep=which(sig_df[,1]==TRUE & sig_df[,2]==FALSE & sig_df[,3]==FALSE & sig_df[,5]==FALSE)
est_comb=est_comb[idx_keep,]


write.csv(est_comb, file="figures/resub/miRglmm_discrepant_miRNA_study89.csv")


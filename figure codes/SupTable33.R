library(ggplot2)

load(file='study46_exact/results/filter_neg1_processed_results.rda')

beta_df=data.frame(results[["beta_hat"]])
beta_df=beta_df[, -c(2,4)]
beta_df=round(beta_df, 2)
CI_LL_df=data.frame(results[["CI_LL"]])[, -c(2,4)]
CI_LL_df=round(CI_LL_df,2)
CI_UL_df=data.frame(results[["CI_UL"]])[, -c(2,4)]
CI_UL_df=round(CI_UL_df,2)
for (ind in seq(1, length(beta_df[1,]))){
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
colnames(est_comb)=colnames(beta_df)
rownames(est_comb)=rownames(beta_df)



p_df=data.frame(results[["pvals"]])
p_df=p_df[, -c(2,4)]
padj_df=matrix(, nrow=dim(p_df)[1], ncol=dim(p_df)[2])
for (ind in seq(1,dim(p_df)[2])){
  padj_df[,ind]=p.adjust(p_df[,ind], method="BH")
}
sig_df=data.frame(padj_df<0.05)
sig_df[padj_df<0.05 & beta_df<0]="sig down"
sig_df[padj_df<0.05 & beta_df>0]="sig up"
colnames(sig_df)=colnames(p_df)
rownames(sig_df)=rownames(data.frame(results[["pvals"]]))

est_comb$miRglmm_p=padj_df[,1]


LRTp=data.frame(results[["LRTp"]])
LRTpadj=data.frame("FDR"=p.adjust(LRTp[,1], method="BH"))
rownames(LRTpadj)=rownames(LRTp)

est_comb=merge(est_comb, LRTpadj, by="row.names")
rownames(est_comb)=est_comb$Row.names
est_comb=est_comb[, -1]

est_comb=est_comb[, c(1,6,7,2:5)]
colnames(est_comb)=c("miRglmm logFC (95% CI)", "fixed effect FDR", "random slope FDR", 
                     "NB GLM logFC (95% CI)", "DESeq2 logFC (95% CI)", "edgeR logFC (95% CI)", 
                     "limma-voom logFC (95% CI)")
all(rownames(est_comb)==rownames(sig_df))

est_comb=est_comb[which(est_comb$`random slope FDR`<1E-7),]


write.csv(est_comb, file="figures/resub/miRglmm_sig_randslope_top10_miRNA.csv")

out=results[["var_comp"]]
out[which(rownames(out) %in% rownames(est_comb)), "random_slope_seq_var"]
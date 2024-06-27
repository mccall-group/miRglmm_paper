library(irr)
library(vcd)
library(ggplot2)

load(file='study 89/results/filter_neg1_processed_results.rda')

beta_df=data.frame(results[["beta_hat"]])
miRNA=rownames(beta_df[beta_df$contrast=="Monocyte vs B lymphocyte",])
miRNA_long=rep(miRNA, 4)
beta_df$miRNA=miRNA_long

p_df=data.frame(results[["pvals"]])
#adjust pvalue for number of miRNA analyzed (dont correct for all contrasts performed)
padj_df=matrix(, nrow=dim(p_df)[1], ncol=dim(p_df)[2]-1)
for (ind in seq(1,dim(p_df)[2]-1)){
  padj_df[which(p_df$contrast=="CD4 T lymphocyte vs B lymphocyte"),ind]=p.adjust(p_df[which(p_df$contrast=="CD4 T lymphocyte vs B lymphocyte"),ind], method="BH")
  padj_df[which(p_df$contrast=="CD8 T lymphocyte vs B lymphocyte"),ind]=p.adjust(p_df[which(p_df$contrast=="CD8 T lymphocyte vs B lymphocyte"),ind], method="BH")
  padj_df[which(p_df$contrast=="Monocyte vs B lymphocyte"),ind]=p.adjust(p_df[which(p_df$contrast=="Monocyte vs B lymphocyte"),ind], method="BH")
  padj_df[which(p_df$contrast=="Natural Killer vs B lymphocyte"),ind]=p.adjust(p_df[which(p_df$contrast=="Natural Killer vs B lymphocyte"),ind], method="BH")
}
sig_df=data.frame(padj_df<0.05)
sig_df$contrast=p_df$contrast
colnames(sig_df)=colnames(data.frame(results[["pvals"]]))
rownames(sig_df)=rownames(data.frame(results[["pvals"]]))
sig_df$miRNA=miRNA_long



#create interaction significant groups of all aggregated methods vs miRglmm
sig_groups_v_miRglmm=data.frame("contrast"=sig_df$contrast, "miRglmnb"=as.character(interaction(sig_df$miRglmm, sig_df$miRglmnb)), "DESeq2"=as.character(interaction(sig_df$miRglmm, sig_df$DESeq2)),
                                "edgeR"=as.character(interaction(sig_df$miRglmm, sig_df$edgeR)), "limmavoom"=as.character(interaction(sig_df$miRglmm, sig_df$limmavoom)))
rownames(sig_groups_v_miRglmm)=rownames(sig_df)



sig_groups_v_miRglmm[sig_groups_v_miRglmm=="TRUE.TRUE"]="significant for both"
sig_groups_v_miRglmm[sig_groups_v_miRglmm=="TRUE.FALSE"]="significant for miRglmm only"
sig_groups_v_miRglmm[sig_groups_v_miRglmm=="FALSE.TRUE"]="significant for aggregation method only"
sig_groups_v_miRglmm[sig_groups_v_miRglmm=="FALSE.FALSE"]="significant for neither"

#combine with the betas to make sure miRNAs align with correct beta
comb_df=transform(merge(beta_df, sig_groups_v_miRglmm, by='row.names'), row.names=Row.names, Row.names=NULL)

idx=which(comb_df$miRglmnb.y=="significant for aggregation method only" & comb_df$DESeq2.y=="significant for aggregation method only" & comb_df$edgeR.y=="significant for aggregation method only" & comb_df$limmavoom.y=="significant for aggregation method only")
table_out=data.frame("miRNA"=comb_df$miRNA[idx], "contrast"=comb_df$contrast.x[idx], "miRglmm logFC"=round(comb_df$miRglmm[idx],2), "NB GLM logFC"=round(comb_df$miRglmnb.x[idx],2),
                     "DESeq2 logFC"=round(comb_df$DESeq2.x[idx],2), "edgeR logFC"= round(comb_df$edgeR.x[idx],2),
                     "limma-voom logFC"=round(comb_df$limmavoom.x[idx],2))

write.csv(table_out, file="figures/SupTable8.csv", row.names=FALSE)
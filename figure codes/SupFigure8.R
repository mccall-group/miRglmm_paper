library(ggplot2)

load(file='study46_exact/results/filter_neg1_processed_results.rda')

beta_df=data.frame(results[["beta_hat"]])
beta_df=beta_df[, -c(2,4)]
est_comb=data.frame("logFC"=beta_df[,1])
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

all(rownames(beta_df)==rownames(p_df))
est_comb$miRglmm_p=padj_df[,1]


LRTp=data.frame(results[["LRTp"]])
LRTpadj=data.frame("FDR"=p.adjust(LRTp[,1], method="BH"))
rownames(LRTpadj)=rownames(LRTp)

est_comb=merge(est_comb, LRTpadj, by="row.names")
rownames(est_comb)=est_comb$Row.names
est_comb=est_comb[, -1]

est_comb$`Differential IsomiR Usage`="Not signficant"
est_comb$`Differential IsomiR Usage`[which(est_comb$FDR<0.05)]="Significant"
est_comb$y=-log10(est_comb$miRglmm_p)
idx=which(est_comb$y==Inf)
est_comb$y[idx]=max(est_comb$y[est_comb$y!=Inf])+0.1

ggplot(est_comb, aes(x=logFC, y=y, color=`Differential IsomiR Usage`))+geom_point()+
  theme_minimal()+scale_color_manual(values=c("black", "red"))+ylim(0,8.5)+
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed')+
  ylab('-log10(FDR)')+xlab('Testes vs Bladder logFC')+
  theme(text=element_text(size=16), 
        legend.text=element_text(size=12))

ggsave("figures/resub/bladder_testes_volcano.tif", plot=last_plot(), device="tiff", width=10, height=7, units="in", dpi=320, bg="white")

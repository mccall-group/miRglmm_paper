library(ggplot2)

load(file='study89_exact/results/filter_neg1_processed_results.rda')
beta_df=data.frame(results[["beta_hat"]])

p_df=data.frame(results[["pvals"]])
padj_df=matrix(, nrow=dim(p_df)[1], ncol=dim(p_df)[2]-1)
for (ind in seq(1,dim(p_df)[2]-1)){
  padj_df[which(p_df$contrast=="CD4 T lymphocyte vs B lymphocyte"),ind]=p.adjust(p_df[which(p_df$contrast=="CD4 T lymphocyte vs B lymphocyte"),ind], method="BH")
  padj_df[which(p_df$contrast=="CD8 T lymphocyte vs B lymphocyte"),ind]=p.adjust(p_df[which(p_df$contrast=="CD8 T lymphocyte vs B lymphocyte"),ind], method="BH")
  padj_df[which(p_df$contrast=="Monocyte vs B lymphocyte"),ind]=p.adjust(p_df[which(p_df$contrast=="Monocyte vs B lymphocyte"),ind], method="BH")
  padj_df[which(p_df$contrast=="Natural Killer vs B lymphocyte"),ind]=p.adjust(p_df[which(p_df$contrast=="Natural Killer vs B lymphocyte"),ind], method="BH")
}

beta_df_all=beta_df
padj_df_all=data.frame(padj_df)
padj_df_all$contrast=p_df[,8]
rownames(padj_df_all)=rownames(p_df)
colnames(padj_df_all)=colnames(p_df)

contrasts=c("Monocyte vs B lymphocyte", "Natural Killer vs B lymphocyte",
            "CD4 T lymphocyte vs B lymphocyte","CD8 T lymphocyte vs B lymphocyte")


p=list()
for (ind in seq(1, length(contrasts))){

contrast_in=contrasts[ind]
padj_df=padj_df_all[which(padj_df_all$contrast==contrast_in), -c(2,4)]
beta_df=beta_df_all[which(beta_df_all$contrast==contrast_in), -c(2,4)]

est_comb=data.frame("logFC"=beta_df[,1])
rownames(est_comb)=rownames(beta_df)

all(rownames(beta_df)==rownames(padj_df))
est_comb$miRglmm_p=padj_df[,1]

if (ind==1){
LRTp=data.frame(results[["LRTp"]])
LRTpadj=data.frame("FDR"=p.adjust(LRTp[,1], method="BH"))
rownames(LRTpadj)=rownames(LRTp)

est_comb=merge(est_comb, LRTpadj, by="row.names")
rownames(est_comb)=est_comb$Row.names
est_comb=est_comb[, -1]
LRTpadj=data.frame("FDR"=est_comb$FDR)
rownames(LRTpadj)=rownames(est_comb)

} else {
  all(rownames(LRTpadj)==rownames(est_comb))
  est_comb$FDR=LRTpadj$FDR
}



est_comb$`Differential IsomiR Usage`="Not signficant"
est_comb$`Differential IsomiR Usage`[which(est_comb$FDR<0.05)]="Significant"
est_comb$y=-log10(est_comb$miRglmm_p)
#
idx=which(est_comb$y>50)
est_comb$y[idx]=50
p[[ind]]=ggplot(est_comb, aes(x=logFC, y=y))+geom_point()+
  theme_minimal()+
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed')+
  ylab('-log10(FDR)')+xlab('logFC')+ggtitle(contrast_in)+
  theme(text=element_text(size=16), 
        legend.text=element_text(size=12))

}
library(ggpubr)
ggarrange(p[[1]], p[[2]], p[[3]], p[[4]], nrow=2, ncol=2)

ggsave("figures/resub/study89_volcano.tif", plot=last_plot(), device="tiff", width=10, height=7, units="in", dpi=320, bg="white")

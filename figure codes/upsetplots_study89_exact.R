library(irr)
library(vcd)
library(ggplot2)
library(tidyverse)

load(file='study89_exact/results/filter_neg1_processed_results.rda')

beta_df=data.frame(results[["beta_hat"]])
p_df=data.frame(results[["pvals"]])
#adjust pvalue for number of miRNA analyzed (dont correct for all contrasts performed)
padj_df=matrix(, nrow=dim(p_df)[1], ncol=dim(p_df)[2]-1)
for (ind in seq(1,dim(p_df)[2]-1)){
  padj_df[which(p_df$contrast=="CD4 T lymphocyte vs B lymphocyte"),ind]=p.adjust(p_df[which(p_df$contrast=="CD4 T lymphocyte vs B lymphocyte"),ind], method="BH")
  padj_df[which(p_df$contrast=="CD8 T lymphocyte vs B lymphocyte"),ind]=p.adjust(p_df[which(p_df$contrast=="CD8 T lymphocyte vs B lymphocyte"),ind], method="BH")
  padj_df[which(p_df$contrast=="Monocyte vs B lymphocyte"),ind]=p.adjust(p_df[which(p_df$contrast=="Monocyte vs B lymphocyte"),ind], method="BH")
  padj_df[which(p_df$contrast=="Natural Killer vs B lymphocyte"),ind]=p.adjust(p_df[which(p_df$contrast=="Natural Killer vs B lymphocyte"),ind], method="BH")
}
beta_df_sub=beta_df[, 1:7]
sig_df=data.frame(padj_df<0.05)
sig_df[padj_df<0.05 & beta_df_sub<0]="sig down"
sig_df[padj_df<0.05 & beta_df_sub>0]="sig up"
sig_df$contrast=p_df$contrast
colnames(sig_df)=colnames(data.frame(results[["pvals"]]))
rownames(sig_df)=rownames(data.frame(results[["pvals"]]))

beta_df_study89=beta_df
results_study89=results
sig_df_study89=sig_df
padj_df_study89=padj_df



sig_df_study89$miRNA=rownames(sig_df_study89)





###########panel B upset plot
library(UpSetR)

#find edgeR significance from fit object
#load(file='bladder_testes_results/filter_neg1_results.rda')
#sig_df_edgeR=data.frame('edgeR'=fits[["edgeR"]][["table"]][["PValue"]]<0.05)
#rownames(sig_df_edgeR)=rownames(fits[["edgeR"]][["table"]])
#sig_df=transform(merge(sig_df, sig_df_edgeR, by='row.names'), row.names=Row.names, Row.names=NULL)

list_Input=list(miRglmm=rownames(sig_df_study89)[sig_df_study89$miRglmm!=FALSE], 
                DESeq2=rownames(sig_df_study89)[sig_df_study89$DESeq2!=FALSE],
                `NB GLM`=rownames(sig_df_study89)[sig_df_study89$miRglmnb!=FALSE],
                `limma-voom`=rownames(sig_df_study89)[sig_df_study89$limmavoom!=FALSE])
myplot=upset(fromList(list_Input), order.by="freq", mainbar.y.label="Number of common significant miRNA",
             sets.x.label="Significant miRNA", text.scale=c(2, 2, 2,2,2,2))

tiff(file="figures/resub/SFX_study89.tif", height=6.74, width=8.9, units="in", res=320)
print(myplot)
dev.off()


list_Input=list(miRglmm=rownames(sig_df_study89)[sig_df_study89$miRglmm!="sig up"], 
                DESeq2=rownames(sig_df_study89)[sig_df_study89$DESeq2!="sig up"],
                `NB GLM`=rownames(sig_df_study89)[sig_df_study89$miRglmnb!="sig up"],
                `limma-voom`=rownames(sig_df_study89)[sig_df_study89$limmavoom!="sig up"])
myplot=upset(fromList(list_Input), order.by="freq", mainbar.y.label="Number of common significant miRNA",
             sets.x.label="Significant miRNA", text.scale=c(2, 2, 2,2,2,2))

tiff(file="figures/resub/SFX_study89_up.tif", height=6.74, width=8.9, units="in", res=320)
print(myplot)
dev.off()

list_Input=list(miRglmm=rownames(sig_df_study89)[sig_df_study89$miRglmm!="sig down"], 
                DESeq2=rownames(sig_df_study89)[sig_df_study89$DESeq2!="sig down"],
                `NB GLM`=rownames(sig_df_study89)[sig_df_study89$miRglmnb!="sig down"],
                `limma-voom`=rownames(sig_df_study89)[sig_df_study89$limmavoom!="sig down"])
myplot=upset(fromList(list_Input), order.by="freq", mainbar.y.label="Number of common significant miRNA",
             sets.x.label="Significant miRNA", text.scale=c(2, 2, 2,2,2,2))

tiff(file="figures/resub/SFX_study89_down.tif", height=6.74, width=8.9, units="in", res=320)
print(myplot)
dev.off()
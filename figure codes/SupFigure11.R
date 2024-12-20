library(irr)
library(vcd)
library(ggplot2)
#library(tidyverse)

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

sig_df_all=sig_df


sig_df=subset(sig_df_all, contrast=="Monocyte vs B lymphocyte")

list_Input=list('miRglmm up'=rownames(sig_df)[sig_df$miRglmm=="sig up"], 'miRglmm down'=rownames(sig_df)[sig_df$miRglmm=="sig down"],
                'DESeq2 up'=rownames(sig_df)[sig_df$DESeq2=="sig up"], 'DESeq2 down'=rownames(sig_df)[sig_df$DESeq2=="sig down"],
                'NB GLM up'=rownames(sig_df)[sig_df$miRglmnb=="sig up"], 'NB GLM down'=rownames(sig_df)[sig_df$miRglmnb=="sig down"],
            #    'edgeR up'=rownames(sig_df)[sig_df$edgeR=="sig up"], 'edgeR down'=rownames(sig_df)[sig_df$edgeR=="sig down"],
                'limma-voom up'=rownames(sig_df)[sig_df$limmavoom=="sig up"], 'limma-voom down'=rownames(sig_df)[sig_df$limmavoom=="sig down"])
myplot=upset(fromList(list_Input), order.by="freq", sets=rev(names(list_Input)), mainbar.y.label="Number of common significant miRNA",
             sets.x.label="Significant miRNA", text.scale=c(2, 2, 2,2,2,2), keep.order=TRUE)
print(myplot)

upset_ind=myplot[["New_data"]]
x1=unlist(list_Input, use.names=TRUE)
x1=x1[!duplicated(x1)]
idx=which((upset_ind$`miRglmm up`==1 | upset_ind$`miRglmm down`==1) &
            upset_ind$`DESeq2 up`==0 & upset_ind$`DESeq2 down`==0 &
            #upset_ind$`edgeR up`==0 & upset_ind$`edgeR down`==0 &
            upset_ind$`limma-voom up`==0 & upset_ind$`limma-voom down`==0 &
            upset_ind$`NB GLM up`==0 & upset_ind$`NB GLM down`==0)
x1[idx]

idx=which(upset_ind$`miRglmm up`==0 & upset_ind$`miRglmm down`==0 &
            upset_ind$`DESeq2 up`==1 & upset_ind$`DESeq2 down`==0 &
            #upset_ind$`edgeR up`==0 & upset_ind$`edgeR down`==0 &
            upset_ind$`limma-voom up`==1 & upset_ind$`limma-voom down`==0 &
            upset_ind$`NB GLM up`==1 & upset_ind$`NB GLM down`==0)
x1[idx]

idx=which(upset_ind$`miRglmm up`==0 & upset_ind$`miRglmm down`==0 &
            upset_ind$`DESeq2 up`==0 & upset_ind$`DESeq2 down`==1 &
            #upset_ind$`edgeR up`==0 & upset_ind$`edgeR down`==0 &
            upset_ind$`limma-voom up`==0 & upset_ind$`limma-voom down`==1 &
            upset_ind$`NB GLM up`==0 & upset_ind$`NB GLM down`==1)
x1[idx]


tiff(file="figures/resub/SFX_study89_monocyte.tif", height=9, width=10, units="in", res=320)
print(myplot)
dev.off()

sig_df=subset(sig_df_all, contrast=="Natural Killer vs B lymphocyte")

list_Input=list('miRglmm up'=rownames(sig_df)[sig_df$miRglmm=="sig up"], 'miRglmm down'=rownames(sig_df)[sig_df$miRglmm=="sig down"],
                'DESeq2 up'=rownames(sig_df)[sig_df$DESeq2=="sig up"], 'DESeq2 down'=rownames(sig_df)[sig_df$DESeq2=="sig down"],
                'NB GLM up'=rownames(sig_df)[sig_df$miRglmnb=="sig up"], 'NB GLM down'=rownames(sig_df)[sig_df$miRglmnb=="sig down"],
                #    'edgeR up'=rownames(sig_df)[sig_df$edgeR=="sig up"], 'edgeR down'=rownames(sig_df)[sig_df$edgeR=="sig down"],
                'limma-voom up'=rownames(sig_df)[sig_df$limmavoom=="sig up"], 'limma-voom down'=rownames(sig_df)[sig_df$limmavoom=="sig down"])
myplot=upset(fromList(list_Input), order.by="freq", sets=rev(names(list_Input)), mainbar.y.label="Number of common significant miRNA",
             sets.x.label="Significant miRNA", text.scale=c(2, 2, 2,2,2,2), keep.order=TRUE)
print(myplot)

upset_ind=myplot[["New_data"]]
x1=unlist(list_Input, use.names=TRUE)
x1=x1[!duplicated(x1)]
idx=which((upset_ind$`miRglmm up`==1 | upset_ind$`miRglmm down`==1) &
            upset_ind$`DESeq2 up`==0 & upset_ind$`DESeq2 down`==0 &
            #upset_ind$`edgeR up`==0 & upset_ind$`edgeR down`==0 &
            upset_ind$`limma-voom up`==0 & upset_ind$`limma-voom down`==0 &
            upset_ind$`NB GLM up`==0 & upset_ind$`NB GLM down`==0)
x1[idx]

idx=which(upset_ind$`miRglmm up`==0 & upset_ind$`miRglmm down`==0 &
            upset_ind$`DESeq2 up`==1 & upset_ind$`DESeq2 down`==0 &
            #upset_ind$`edgeR up`==0 & upset_ind$`edgeR down`==0 &
            upset_ind$`limma-voom up`==1 & upset_ind$`limma-voom down`==0 &
            upset_ind$`NB GLM up`==1 & upset_ind$`NB GLM down`==0)
x1[idx]

idx=which(upset_ind$`miRglmm up`==0 & upset_ind$`miRglmm down`==0 &
            upset_ind$`DESeq2 up`==0 & upset_ind$`DESeq2 down`==1 &
            #upset_ind$`edgeR up`==0 & upset_ind$`edgeR down`==0 &
            upset_ind$`limma-voom up`==0 & upset_ind$`limma-voom down`==1 &
            upset_ind$`NB GLM up`==0 & upset_ind$`NB GLM down`==1)
x1[idx]



tiff(file="figures/resub/SFX_study89_NK.tif", height=9, width=10, units="in", res=320)
print(myplot)
dev.off()

sig_df=subset(sig_df_all, contrast=="CD4 T lymphocyte vs B lymphocyte")

list_Input=list('miRglmm up'=rownames(sig_df)[sig_df$miRglmm=="sig up"], 'miRglmm down'=rownames(sig_df)[sig_df$miRglmm=="sig down"],
                'DESeq2 up'=rownames(sig_df)[sig_df$DESeq2=="sig up"], 'DESeq2 down'=rownames(sig_df)[sig_df$DESeq2=="sig down"],
                'NB GLM up'=rownames(sig_df)[sig_df$miRglmnb=="sig up"], 'NB GLM down'=rownames(sig_df)[sig_df$miRglmnb=="sig down"],
                #    'edgeR up'=rownames(sig_df)[sig_df$edgeR=="sig up"], 'edgeR down'=rownames(sig_df)[sig_df$edgeR=="sig down"],
                'limma-voom up'=rownames(sig_df)[sig_df$limmavoom=="sig up"], 'limma-voom down'=rownames(sig_df)[sig_df$limmavoom=="sig down"])
myplot=upset(fromList(list_Input), order.by="freq", sets=rev(names(list_Input)), mainbar.y.label="Number of common significant miRNA",
             sets.x.label="Significant miRNA", text.scale=c(2, 2, 2,2,2,2), keep.order=TRUE)
print(myplot)

upset_ind=myplot[["New_data"]]
x1=unlist(list_Input, use.names=TRUE)
x1=x1[!duplicated(x1)]
idx=which((upset_ind$`miRglmm up`==1 | upset_ind$`miRglmm down`==1) &
            upset_ind$`DESeq2 up`==0 & upset_ind$`DESeq2 down`==0 &
            #upset_ind$`edgeR up`==0 & upset_ind$`edgeR down`==0 &
            upset_ind$`limma-voom up`==0 & upset_ind$`limma-voom down`==0 &
            upset_ind$`NB GLM up`==0 & upset_ind$`NB GLM down`==0)
x1[idx]

idx=which(upset_ind$`miRglmm up`==0 & upset_ind$`miRglmm down`==0 &
            upset_ind$`DESeq2 up`==1 & upset_ind$`DESeq2 down`==0 &
            #upset_ind$`edgeR up`==0 & upset_ind$`edgeR down`==0 &
            upset_ind$`limma-voom up`==1 & upset_ind$`limma-voom down`==0 &
            upset_ind$`NB GLM up`==1 & upset_ind$`NB GLM down`==0)
x1[idx]

idx=which(upset_ind$`miRglmm up`==0 & upset_ind$`miRglmm down`==0 &
            upset_ind$`DESeq2 up`==0 & upset_ind$`DESeq2 down`==1 &
            #upset_ind$`edgeR up`==0 & upset_ind$`edgeR down`==0 &
            upset_ind$`limma-voom up`==0 & upset_ind$`limma-voom down`==1 &
            upset_ind$`NB GLM up`==0 & upset_ind$`NB GLM down`==1)
x1[idx]



tiff(file="figures/resub/SFX_study89_CD4.tif", height=9, width=10, units="in", res=320)
print(myplot)
dev.off()

sig_df=subset(sig_df_all, contrast=="CD8 T lymphocyte vs B lymphocyte")

list_Input=list('miRglmm up'=rownames(sig_df)[sig_df$miRglmm=="sig up"], 'miRglmm down'=rownames(sig_df)[sig_df$miRglmm=="sig down"],
                'DESeq2 up'=rownames(sig_df)[sig_df$DESeq2=="sig up"], 'DESeq2 down'=rownames(sig_df)[sig_df$DESeq2=="sig down"],
                'NB GLM up'=rownames(sig_df)[sig_df$miRglmnb=="sig up"], 'NB GLM down'=rownames(sig_df)[sig_df$miRglmnb=="sig down"],
                #    'edgeR up'=rownames(sig_df)[sig_df$edgeR=="sig up"], 'edgeR down'=rownames(sig_df)[sig_df$edgeR=="sig down"],
                'limma-voom up'=rownames(sig_df)[sig_df$limmavoom=="sig up"], 'limma-voom down'=rownames(sig_df)[sig_df$limmavoom=="sig down"])
myplot=upset(fromList(list_Input), order.by="freq", sets=rev(names(list_Input)), mainbar.y.label="Number of common significant miRNA",
             sets.x.label="Significant miRNA", text.scale=c(2, 2, 2,2,2,2), keep.order=TRUE)
print(myplot)

upset_ind=myplot[["New_data"]]
x1=unlist(list_Input, use.names=TRUE)
x1=x1[!duplicated(x1)]
idx=which((upset_ind$`miRglmm up`==1 | upset_ind$`miRglmm down`==1) &
            upset_ind$`DESeq2 up`==0 & upset_ind$`DESeq2 down`==0 &
            #upset_ind$`edgeR up`==0 & upset_ind$`edgeR down`==0 &
            upset_ind$`limma-voom up`==0 & upset_ind$`limma-voom down`==0 &
            upset_ind$`NB GLM up`==0 & upset_ind$`NB GLM down`==0)
x1[idx]

idx=which(upset_ind$`miRglmm up`==0 & upset_ind$`miRglmm down`==0 &
            upset_ind$`DESeq2 up`==1 & upset_ind$`DESeq2 down`==0 &
            #upset_ind$`edgeR up`==0 & upset_ind$`edgeR down`==0 &
            upset_ind$`limma-voom up`==1 & upset_ind$`limma-voom down`==0 &
            upset_ind$`NB GLM up`==1 & upset_ind$`NB GLM down`==0)
x1[idx]

idx=which(upset_ind$`miRglmm up`==0 & upset_ind$`miRglmm down`==0 &
            upset_ind$`DESeq2 up`==0 & upset_ind$`DESeq2 down`==1 &
            #upset_ind$`edgeR up`==0 & upset_ind$`edgeR down`==0 &
            upset_ind$`limma-voom up`==0 & upset_ind$`limma-voom down`==1 &
            upset_ind$`NB GLM up`==0 & upset_ind$`NB GLM down`==1)
x1[idx]



tiff(file="figures/resub/SFX_study89_CD8.tif", height=9, width=10, units="in", res=320)
print(myplot)
dev.off()


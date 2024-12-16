library(irr)
library(vcd)
library(ggplot2)

load(file='study46_exact/results/filter_neg1_processed_results.rda')

beta_df=data.frame(results[["beta_hat"]])
p_df=data.frame(results[["pvals"]])
padj_df=matrix(, nrow=dim(p_df)[1], ncol=dim(p_df)[2])
for (ind in seq(1,dim(p_df)[2])){
  padj_df[,ind]=p.adjust(p_df[,ind], method="BH")
}
sig_df=data.frame(padj_df<0.05)
sig_df[padj_df<0.05 & beta_df<0]="sig down"
sig_df[padj_df<0.05 & beta_df>0]="sig up"
colnames(sig_df)=colnames(data.frame(results[["pvals"]]))
rownames(sig_df)=rownames(data.frame(results[["pvals"]]))


beta_df_study46=beta_df
beta_df_study46$miRNA=rownames(beta_df_study46)
results_study46=results
sig_df_study46=sig_df
padj_df_study46=padj_df


#load miRblood results to compare
load(file='GTex_bladder_testes_exact/results/nofilter_processed_results.rda')

beta_df=data.frame(results[["beta_hat"]])
p_df=data.frame(results[["pvals"]])
padj_df=matrix(, nrow=dim(p_df)[1], ncol=dim(p_df)[2])
for (ind in seq(1,dim(p_df)[2])){
  padj_df[,ind]=p.adjust(p_df[,ind], method="BH")
}
sig_df=data.frame(padj_df<0.05)
sig_df[padj_df<0.05 & beta_df<0]="sig down"
sig_df[padj_df<0.05 & beta_df>0]="sig up"
colnames(sig_df)=colnames(data.frame(results[["pvals"]]))
rownames(sig_df)=rownames(data.frame(results[["pvals"]]))


beta_df_gtex=beta_df
beta_df_gtex$miRNA=rownames(beta_df_gtex)
results_gtex=results
sig_df_gtex=sig_df
padj_df_gtex=padj_df

sig_df_gtex$miRNA=rownames(sig_df_gtex)
sig_df_study46$miRNA=rownames(sig_df_study46)


miRglmm_study46=beta_df_study46[, c("miRglmm" ,"miRNA")]
miRglmm_gtex=beta_df_gtex[, c("miRglmm" ,"miRNA")]

miRglmm_merge=merge(miRglmm_study46, miRglmm_gtex, by=c("miRNA"))
miRglmm_correlation=cor(miRglmm_merge$miRglmm.x, miRglmm_merge$miRglmm.y, method="spearman")
pearson_miRglmm_correlation=cor(miRglmm_merge$miRglmm.x, miRglmm_merge$miRglmm.y, method="pearson")
p1=ggplot(miRglmm_merge, aes(x=miRglmm.x, y=miRglmm.y))+geom_point()+
  theme_minimal()+
  xlab('Li logFC estimate')+ylab('GTEx logFC estimate')+ggtitle('miRglmm')+
  theme(plot.title=element_text(size=12, hjust=0.5))

miRglmnb_study46=beta_df_study46[, c("miRglmnb" ,"miRNA")]
miRglmnb_gtex=beta_df_gtex[, c("miRglmnb" ,"miRNA")]

miRglmnb_merge=merge(miRglmnb_study46, miRglmnb_gtex, by=c("miRNA"))
miRglmnb_correlation=cor(miRglmnb_merge$miRglmnb.x, miRglmnb_merge$miRglmnb.y, method="spearman")
pearson_miRglmnb_correlation=cor(miRglmnb_merge$miRglmnb.x, miRglmnb_merge$miRglmnb.y, method="pearson")
p2=ggplot(miRglmnb_merge, aes(x=miRglmnb.x, y=miRglmnb.y))+geom_point()+
  theme_minimal()+
  xlab('Li logFC estimate')+ylab('GTEx logFC estimate')+ggtitle('NB GLM')+
  theme(plot.title=element_text(size=12, hjust=0.5))

deseq_study46=beta_df_study46[, c("DESeq2" ,"miRNA")]
deseq_gtex=beta_df_gtex[, c("DESeq2" ,"miRNA")]

deseq_merge=merge(deseq_study46, deseq_gtex, by=c("miRNA"))
deseq_correlation=cor(deseq_merge$DESeq2.x, deseq_merge$DESeq2.y, method="spearman")
pearson_deseq_correlation=cor(deseq_merge$DESeq2.x, deseq_merge$DESeq2.y, method="pearson")
p3=ggplot(deseq_merge, aes(x=DESeq2.x, y=DESeq2.y))+geom_point()+
  theme_minimal()+
  xlab('Li logFC estimate')+ylab('GTEx logFC estimate')+ggtitle('DESeq2')+
  theme(plot.title=element_text(size=12, hjust=0.5))

edgeR_study46=beta_df_study46[, c("edgeR" ,"miRNA")]
edgeR_gtex=beta_df_gtex[, c("edgeR" ,"miRNA")]

edgeR_merge=merge(edgeR_study46, edgeR_gtex, by=c("miRNA"))
edgeR_correlation=cor(edgeR_merge$edgeR.x, edgeR_merge$edgeR.y, method="spearman")
pearson_edgeR_correlation=cor(edgeR_merge$edgeR.x, edgeR_merge$edgeR.y, method="pearson")
p4=ggplot(edgeR_merge, aes(x=edgeR.x, y=edgeR.y))+geom_point()+
  theme_minimal()+
  xlab('Li logFC estimate')+ylab('GTEx logFC estimate')+ggtitle('edgeR')+
  theme(plot.title=element_text(size=12, hjust=0.5))

limmavoom_study46=beta_df_study46[, c("limmavoom" ,"miRNA")]
limmavoom_gtex=beta_df_gtex[, c("limmavoom" ,"miRNA")]

limmavoom_merge=merge(limmavoom_study46, limmavoom_gtex, by=c("miRNA"))
limmavoom_correlation=cor(limmavoom_merge$limmavoom.x, limmavoom_merge$limmavoom.y, method="spearman")
pearson_limmavoom_correlation=cor(limmavoom_merge$limmavoom.x, limmavoom_merge$limmavoom.y, method="pearson")
p5=ggplot(limmavoom_merge, aes(x=limmavoom.x, y=limmavoom.y))+geom_point()+
  theme_minimal()+
  xlab('Li logFC estimate')+ylab('GTEx logFC estimate')+ggtitle('limma-voom')+
  theme(plot.title=element_text(size=12, hjust=0.5))

library(ggpubr)
ggarrange(p1, p2, p3, p4, p5, ncol=2, nrow=3)

corr_all=data.frame("spearman"=rbind(miRglmm_correlation, miRglmnb_correlation, deseq_correlation, edgeR_correlation, limmavoom_correlation))
rownames(corr_all)=c("miRglmm", "NB GLM", "DESeq2", "edgeR", "limma-voom")

pearson_corr_all=data.frame("pearson"=rbind(pearson_miRglmm_correlation, pearson_miRglmnb_correlation, pearson_deseq_correlation, pearson_edgeR_correlation, pearson_limmavoom_correlation))
corr_all$pearson=pearson_corr_all[,1]

sigdf_miRglmm_merge=merge(sig_df_study46[, c("miRglmm" ,"miRNA")], sig_df_gtex[, c("miRglmm" ,"miRNA")], by=c("miRNA"))
sigdf_miRglmnb_merge=merge(sig_df_study46[, c("miRglmnb" ,"miRNA")], sig_df_gtex[, c("miRglmnb" ,"miRNA")], by=c("miRNA"))
sigdf_deseq_merge=merge(sig_df_study46[, c("DESeq2" ,"miRNA")], sig_df_gtex[, c("DESeq2" ,"miRNA")], by=c("miRNA"))
sigdf_edgeR_merge=merge(sig_df_study46[, c("edgeR" ,"miRNA")], sig_df_gtex[, c("edgeR" ,"miRNA")], by=c("miRNA"))
sigdf_limmavoom_merge=merge(sig_df_study46[, c("limmavoom" ,"miRNA")], sig_df_gtex[, c("limmavoom" ,"miRNA")], by=c("miRNA"))


# miRglmm_kappa=Kappa(table(sigdf_miRglmm_merge$miRglmm.x, sigdf_miRglmm_merge$miRglmm.y))[["Unweighted"]][["value"]]
# miRglmnb_kappa=Kappa(table(sigdf_miRglmnb_merge$miRglmnb.x, sigdf_miRglmnb_merge$miRglmnb.y))[["Unweighted"]][["value"]]
# deseq_kappa=Kappa(table(sigdf_deseq_merge$DESeq2.x, sigdf_deseq_merge$DESeq2.y))[["Unweighted"]][["value"]]
# edgeR_kappa=Kappa(table(sigdf_edgeR_merge$edgeR.x, sigdf_edgeR_merge$edgeR.y))[["Unweighted"]][["value"]]
# limmavoom_kappa=Kappa(table(sigdf_limmavoom_merge$limmavoom.x, sigdf_limmavoom_merge$limmavoom.y))[["Unweighted"]][["value"]]
# 
# 
# 
# kappa_all=data.frame("kappa"=rbind(miRglmm_kappa, miRglmnb_kappa, deseq_kappa, edgeR_kappa, limmavoom_kappa))
# corr_all$kappa=kappa_all[,1]

sig_agree=data.frame(matrix(data=0, nrow=dim(sig_df)[1]))
rownames(sig_agree)=rownames(sig_df)
idx=which(sigdf_miRglmm_merge$miRglmm.x=="sig up" & sigdf_miRglmm_merge$miRglmm.y=="sig up")
sig_agree[idx,1]=1
idx=which(sigdf_miRglmm_merge$miRglmm.x=="sig down" & sigdf_miRglmm_merge$miRglmm.y=="sig down")
sig_agree[idx,1]=1
colnames(sig_agree)="miRglmm"
sig_agree$miRglmnb=0
idx=which(sigdf_miRglmnb_merge$miRglmnb.x=="sig up" & sigdf_miRglmnb_merge$miRglmnb.y=="sig up")
sig_agree$miRglmnb[idx]=1
idx=which(sigdf_miRglmnb_merge$miRglmnb.x=="sig down" & sigdf_miRglmnb_merge$miRglmnb.y=="sig down")
sig_agree$miRglmnb[idx]=1
sig_agree$DESeq2=0
idx=which(sigdf_deseq_merge$DESeq2.x=="sig up" & sigdf_deseq_merge$DESeq2.y=="sig up")
sig_agree$DESeq2[idx]=1
idx=which(sigdf_deseq_merge$DESeq2.x=="sig down" & sigdf_deseq_merge$DESeq2.y=="sig down")
sig_agree$DESeq2[idx]=1
sig_agree$edgeR=0
idx=which(sigdf_edgeR_merge$edgeR.x=="sig up" & sigdf_edgeR_merge$edgeR.y=="sig up")
sig_agree$edgeR[idx]=1
idx=which(sigdf_edgeR_merge$edgeR.x=="sig down" & sigdf_edgeR_merge$edgeR.y=="sig down")
sig_agree$edgeR[idx]=1
sig_agree$limmavoom=0
idx=which(sigdf_limmavoom_merge$limmavoom.x=="sig up" & sigdf_limmavoom_merge$limmavoom.y=="sig up")
sig_agree$limmavoom[idx]=1
idx=which(sigdf_limmavoom_merge$limmavoom.x=="sig down" & sigdf_limmavoom_merge$limmavoom.y=="sig down")
sig_agree$limmavoom[idx]=1

sig_df=sig_agree

###########panel B upset plot
library(UpSetR)

list_Input=list(miRglmm=rownames(sig_df)[sig_df$miRglmm==1], DESeq2=rownames(sig_df)[sig_df$DESeq2==1],`NB GLM`=rownames(sig_df)[sig_df$miRglmnb==1], 
                edgeR=rownames(sig_df)[sig_df$edgeR==1], `limma-voom`=rownames(sig_df)[sig_df$limmavoom==1])
myplot=upset(fromList(list_Input), order.by="freq", mainbar.y.label="Number of common significant miRNA",
             sets.x.label="Significant miRNA", text.scale=c(2, 2, 2,2,2,2))

tiff(file="figures/resub/figure3B.tif", height=6.74, width=8.9, units="in", res=320)
print(myplot)
dev.off()

list_Input=list(miRglmm=rownames(sig_df_study46)[sig_df_study46$miRglmm!=FALSE], DESeq2=rownames(sig_df_study46)[sig_df_study46$DESeq2!=FALSE],
                `NB GLM`=rownames(sig_df_study46)[sig_df_study46$miRglmnb!=FALSE],
                edgeR=rownames(sig_df_study46)[sig_df_study46$edgeR!=FALSE], `limma-voom`=rownames(sig_df_study46)[sig_df_study46$limmavoom!=FALSE])
myplot=upset(fromList(list_Input), order.by="freq", mainbar.y.label="Number of common significant miRNA",
             sets.x.label="Significant miRNA", text.scale=c(2, 2, 2,2,2,2))

tiff(file="figures/resub/figure3B_study46.tif", height=6.74, width=8.9, units="in", res=320)
print(myplot)
dev.off()

list_Input=list(miRglmm=rownames(sig_df_gtex)[sig_df_gtex$miRglmm!=FALSE], DESeq2=rownames(sig_df_gtex)[sig_df_gtex$DESeq2!=FALSE],
                `NB GLM`=rownames(sig_df_gtex)[sig_df_gtex$miRglmnb!=FALSE],
                edgeR=rownames(sig_df_gtex)[sig_df_gtex$edgeR!=FALSE], `limma-voom`=rownames(sig_df_gtex)[sig_df_gtex$limmavoom!=FALSE])
myplot=upset(fromList(list_Input), order.by="freq", mainbar.y.label="Number of common significant miRNA",
             sets.x.label="Significant miRNA", text.scale=c(2, 2, 2,2,2,2))
tiff(file="figures/resub/figure3B_gtex.tif", height=6.74, width=8.9, units="in", res=320)
print(myplot)
dev.off()

list_Input=list(miRglmm=rownames(sig_df_gtex)[sig_df_gtex$miRglmm=="sig up"], DESeq2=rownames(sig_df_gtex)[sig_df_gtex$DESeq2=="sig up"],
                `NB GLM`=rownames(sig_df_gtex)[sig_df_gtex$miRglmnb=="sig up"],
                edgeR=rownames(sig_df_gtex)[sig_df_gtex$edgeR=="sig up"], `limma-voom`=rownames(sig_df_gtex)[sig_df_gtex$limmavoom=="sig up"])
myplot=upset(fromList(list_Input), order.by="freq", mainbar.y.label="Number of common significant miRNA",
             sets.x.label="Significant miRNA", text.scale=c(2, 2, 2,2,2,2))
tiff(file="figures/resub/figure3B_gtex_up.tif", height=6.74, width=8.9, units="in", res=320)
print(myplot)
dev.off()

list_Input=list(miRglmm=rownames(sig_df_gtex)[sig_df_gtex$miRglmm=="sig down"], DESeq2=rownames(sig_df_gtex)[sig_df_gtex$DESeq2=="sig down"],
                `NB GLM`=rownames(sig_df_gtex)[sig_df_gtex$miRglmnb=="sig down"],
                edgeR=rownames(sig_df_gtex)[sig_df_gtex$edgeR=="sig down"], `limma-voom`=rownames(sig_df_gtex)[sig_df_gtex$limmavoom=="sig down"])
myplot=upset(fromList(list_Input), order.by="freq", mainbar.y.label="Number of common significant miRNA",
             sets.x.label="Significant miRNA", text.scale=c(2, 2, 2,2,2,2))
tiff(file="figures/resub/figure3B_gtex_down.tif", height=6.74, width=8.9, units="in", res=320)
print(myplot)
dev.off()
# upset_ind=myplot[["New_data"]]
# x1=unlist(list_Input, use.names=TRUE)
# x1=x1[!duplicated(x1)]
# idx=which(upset_ind$miRglmm==1 & upset_ind$DESeq2==0 & upset_ind$edgeR==0 & upset_ind$`limma-voom`==0 & upset_ind$`NB GLM`==0)
# x1[idx]
# idx=which(upset_ind$miRglmm==0 & upset_ind$DESeq2==1 & upset_ind$edgeR==1 & upset_ind$`limma-voom`==1 & upset_ind$`NB GLM`==1)
# x1[idx]


########## panel C 
load(file='study46_exact/results/filter_neg1_processed_results.rda')
LRTp=data.frame(results[["LRTp"]])
p6=ggplot(LRTp, aes(x=LRTp))+geom_histogram(color="black", fill="gray", bins=50)+xlab('Likelihood Ratio Test p-value')+ylab('number of miRNA')+
  theme(axis.title.x=element_text(size=15), axis.title.y=element_text(size=15),
        axis.text.x=element_text(size=15), axis.text.y=element_text(size=15))
#ggsave("figures/figure4C.tif", plot=last_plot(), device="tiff", scale=6, width=40, height=20, units="mm", dpi=320, bg="white")


library(ggpubr)
ggarrange(ggarrange(p1, p2, p3, p4, p5, ncol=3, nrow=2), ggarrange(NULL,p6, nrow=2), ncol=2)
ggsave("figures/resub/figure3.tif", plot=last_plot(), device="tiff", width=12.8, height=6.64, units="in", dpi=320, bg="white")

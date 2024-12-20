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


sigdf_miRglmm_merge=merge(sig_df_study46[, c("miRglmm" ,"miRNA")], sig_df_gtex[, c("miRglmm" ,"miRNA")], by=c("miRNA"))
sigdf_miRglmnb_merge=merge(sig_df_study46[, c("miRglmnb" ,"miRNA")], sig_df_gtex[, c("miRglmnb" ,"miRNA")], by=c("miRNA"))
sigdf_deseq_merge=merge(sig_df_study46[, c("DESeq2" ,"miRNA")], sig_df_gtex[, c("DESeq2" ,"miRNA")], by=c("miRNA"))
sigdf_edgeR_merge=merge(sig_df_study46[, c("edgeR" ,"miRNA")], sig_df_gtex[, c("edgeR" ,"miRNA")], by=c("miRNA"))
sigdf_limmavoom_merge=merge(sig_df_study46[, c("limmavoom" ,"miRNA")], sig_df_gtex[, c("limmavoom" ,"miRNA")], by=c("miRNA"))


miRglmm_merge=merge(miRglmm_study46, miRglmm_gtex, by=c("miRNA"))
miRglmm_merge=merge(miRglmm_merge, sigdf_miRglmm_merge, by=c("miRNA"))
colnames(miRglmm_merge)=c("miRNA", "study46_beta", "GTEx_beta", "sig_study46", "sig_GTEx")
miRglmm_merge$sig_agree="Not significant in one or both datasets"
miRglmm_merge$sig_agree[which((miRglmm_merge$sig_study46=="sig up" & miRglmm_merge$sig_GTEx=="sig up") | 
                                (miRglmm_merge$sig_study46=="sig down" & miRglmm_merge$sig_GTEx=="sig down"))]="Significant in both datasets"
miRglmm_correlation=cor(miRglmm_merge$study46_beta, miRglmm_merge$GTEx_beta, method="spearman")
pearson_miRglmm_correlation=cor(miRglmm_merge$study46_beta, miRglmm_merge$GTEx_beta, method="pearson")
p1=ggplot(miRglmm_merge, aes(x=study46_beta, y=GTEx_beta, color=sig_agree))+geom_point()+
  theme_minimal()+scale_color_manual(values=c("black", "red"))+
  annotate("text", x=Inf, y=-Inf, label=paste("R =", format(round(pearson_miRglmm_correlation,2), nsmall=2)), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  xlab('Li logFC estimate')+ylab('GTEx logFC estimate')+ggtitle('miRglmm')+
  theme(plot.title=element_text(size=20, hjust=0.5), 
        text=element_text(size=16), 
        legend.text=element_text(size=12), 
        legend.title=element_blank(), 
        legend.position="none")

p1_legend=ggplot(miRglmm_merge, aes(x=study46_beta, y=GTEx_beta, color=sig_agree))+geom_point()+
  theme_minimal()+scale_color_manual(values=c("black", "red"))+
  annotate("text", x=Inf, y=-Inf, label=paste("R =", format(round(pearson_miRglmm_correlation,2), nsmall=2)), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  xlab('Li logFC estimate')+ylab('GTEx logFC estimate')+ggtitle('miRglmm')+
  theme(plot.title=element_text(size=20, hjust=0.5), 
        text=element_text(size=16), 
        legend.text=element_text(size=12), 
        legend.title=element_blank())
ggsave("figures/resub/figure3_legend.tif", plot=last_plot(), device="tiff", width=12.8, height=6.64, units="in", dpi=320, bg="white")


miRglmnb_study46=beta_df_study46[, c("miRglmnb" ,"miRNA")]
miRglmnb_gtex=beta_df_gtex[, c("miRglmnb" ,"miRNA")]

miRglmnb_merge=merge(miRglmnb_study46, miRglmnb_gtex, by=c("miRNA"))
miRglmnb_merge=merge(miRglmnb_merge, sigdf_miRglmnb_merge, by=c("miRNA"))
colnames(miRglmnb_merge)=c("miRNA", "study46_beta", "GTEx_beta", "sig_study46", "sig_GTEx")
miRglmnb_merge$sig_agree="Not significant in one or both datasets"
miRglmnb_merge$sig_agree[which((miRglmnb_merge$sig_study46=="sig up" & miRglmnb_merge$sig_GTEx=="sig up") | 
                                (miRglmnb_merge$sig_study46=="sig down" & miRglmnb_merge$sig_GTEx=="sig down"))]="Significant in both datasets"
miRglmnb_correlation=cor(miRglmnb_merge$study46_beta, miRglmnb_merge$GTEx_beta, method="spearman")
pearson_miRglmnb_correlation=cor(miRglmnb_merge$study46_beta, miRglmnb_merge$GTEx_beta, method="pearson")
p2=ggplot(miRglmnb_merge, aes(x=study46_beta, y=GTEx_beta, color=sig_agree))+geom_point()+
  theme_minimal()+scale_color_manual(values=c("black", "red"))+
  annotate("text", x=Inf, y=-Inf, label=paste("R =", format(round(pearson_miRglmnb_correlation,2), nsmall=2)), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  xlab('Li logFC estimate')+ylab('GTEx logFC estimate')+ggtitle('NB GLM')+
  theme(plot.title=element_text(size=20, hjust=0.5), 
        text=element_text(size=16), 
        legend.text=element_text(size=12), 
        legend.title=element_blank(), 
        legend.position="none")

deseq_study46=beta_df_study46[, c("DESeq2" ,"miRNA")]
deseq_gtex=beta_df_gtex[, c("DESeq2" ,"miRNA")]

deseq_merge=merge(deseq_study46, deseq_gtex, by=c("miRNA"))
deseq_merge=merge(deseq_merge, sigdf_deseq_merge, by=c("miRNA"))
colnames(deseq_merge)=c("miRNA", "study46_beta", "GTEx_beta", "sig_study46", "sig_GTEx")
deseq_merge$sig_agree="Not significant in one or both datasets"
deseq_merge$sig_agree[which((deseq_merge$sig_study46=="sig up" & deseq_merge$sig_GTEx=="sig up") | 
                                 (deseq_merge$sig_study46=="sig down" & deseq_merge$sig_GTEx=="sig down"))]="Significant in both datasets"
deseq_correlation=cor(deseq_merge$study46_beta, deseq_merge$GTEx_beta, method="spearman")
pearson_deseq_correlation=cor(deseq_merge$study46_beta, deseq_merge$GTEx_beta, method="pearson")
p3=ggplot(deseq_merge, aes(x=study46_beta, y=GTEx_beta, color=sig_agree))+geom_point()+
  theme_minimal()+scale_color_manual(values=c("black", "red"))+
  annotate("text", x=Inf, y=-Inf, label=paste("R =", format(round(pearson_deseq_correlation,2), nsmall=2)), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  xlab('Li logFC estimate')+ylab('GTEx logFC estimate')+ggtitle('DESeq2')+
  theme(plot.title=element_text(size=20, hjust=0.5), 
        text=element_text(size=16), 
        legend.text=element_text(size=12), 
        legend.title=element_blank(), 
        legend.position="none")

edgeR_study46=beta_df_study46[, c("edgeR" ,"miRNA")]
edgeR_gtex=beta_df_gtex[, c("edgeR" ,"miRNA")]

edgeR_merge=merge(edgeR_study46, edgeR_gtex, by=c("miRNA"))
edgeR_merge=merge(edgeR_merge, sigdf_edgeR_merge, by=c("miRNA"))
colnames(edgeR_merge)=c("miRNA", "study46_beta", "GTEx_beta", "sig_study46", "sig_GTEx")
edgeR_merge$sig_agree="Not significant in one or both datasets"
edgeR_merge$sig_agree[which((edgeR_merge$sig_study46=="sig up" & edgeR_merge$sig_GTEx=="sig up") | 
                              (edgeR_merge$sig_study46=="sig down" & edgeR_merge$sig_GTEx=="sig down"))]="Significant in both datasets"
edgeR_correlation=cor(edgeR_merge$study46_beta, edgeR_merge$GTEx_beta, method="spearman")
pearson_edgeR_correlation=cor(edgeR_merge$study46_beta, edgeR_merge$GTEx_beta, method="pearson")
p4=ggplot(edgeR_merge, aes(x=study46_beta, y=GTEx_beta, color=sig_agree))+geom_point()+
  theme_minimal()+scale_color_manual(values=c("black", "red"))+
  annotate("text", x=Inf, y=-Inf, label=paste("R =", format(round(pearson_edgeR_correlation,2), nsmall=2)), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  xlab('Li logFC estimate')+ylab('GTEx logFC estimate')+ggtitle('edgeR')+
  theme(plot.title=element_text(size=20, hjust=0.5), 
        text=element_text(size=16), 
        legend.text=element_text(size=12), 
        legend.title=element_blank(), 
        legend.position="none")

limmavoom_study46=beta_df_study46[, c("limmavoom" ,"miRNA")]
limmavoom_gtex=beta_df_gtex[, c("limmavoom" ,"miRNA")]

limmavoom_merge=merge(limmavoom_study46, limmavoom_gtex, by=c("miRNA"))
limmavoom_merge=merge(limmavoom_merge, sigdf_limmavoom_merge, by=c("miRNA"))
colnames(limmavoom_merge)=c("miRNA", "study46_beta", "GTEx_beta", "sig_study46", "sig_GTEx")
limmavoom_merge$sig_agree="Not significant in one or both datasets"
limmavoom_merge$sig_agree[which((limmavoom_merge$sig_study46=="sig up" & limmavoom_merge$sig_GTEx=="sig up") | 
                              (limmavoom_merge$sig_study46=="sig down" & limmavoom_merge$sig_GTEx=="sig down"))]="Significant in both datasets"
limmavoom_correlation=cor(limmavoom_merge$study46_beta, limmavoom_merge$GTEx_beta, method="spearman")
pearson_limmavoom_correlation=cor(limmavoom_merge$study46_beta, limmavoom_merge$GTEx_beta, method="pearson")
p5=ggplot(limmavoom_merge, aes(x=study46_beta, y=GTEx_beta, color=sig_agree))+geom_point()+
  theme_minimal()+scale_color_manual(values=c("black", "red"))+
  annotate("text", x=Inf, y=-Inf, label=paste("R =", format(round(pearson_limmavoom_correlation,2), nsmall=2)), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  xlab('Li logFC estimate')+ylab('GTEx logFC estimate')+ggtitle('limma-voom')+
  theme(plot.title=element_text(size=20, hjust=0.5), 
        text=element_text(size=16), 
        legend.text=element_text(size=12), 
        legend.title=element_blank(), 
        legend.position="none")

library(ggpubr)
ggarrange(p1, p2, p3, p4, p5, ncol=2, nrow=3)
ggsave("figures/resub/study46_vs_GTEx_scatteronly.tif", plot=last_plot(), device="tiff", width=10, height=12, units="in", dpi=320, bg="white")


corr_all=data.frame("spearman"=rbind(miRglmm_correlation, miRglmnb_correlation, deseq_correlation, edgeR_correlation, limmavoom_correlation))
rownames(corr_all)=c("miRglmm", "NB GLM", "DESeq2", "edgeR", "limma-voom")

pearson_corr_all=data.frame("pearson"=rbind(pearson_miRglmm_correlation, pearson_miRglmnb_correlation, pearson_deseq_correlation, pearson_edgeR_correlation, pearson_limmavoom_correlation))
corr_all$pearson=pearson_corr_all[,1]



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
sig_agree[idx,1]=2
colnames(sig_agree)="miRglmm"
sig_agree$miRglmnb=0
idx=which(sigdf_miRglmnb_merge$miRglmnb.x=="sig up" & sigdf_miRglmnb_merge$miRglmnb.y=="sig up")
sig_agree$miRglmnb[idx]=1
idx=which(sigdf_miRglmnb_merge$miRglmnb.x=="sig down" & sigdf_miRglmnb_merge$miRglmnb.y=="sig down")
sig_agree$miRglmnb[idx]=2
sig_agree$DESeq2=0
idx=which(sigdf_deseq_merge$DESeq2.x=="sig up" & sigdf_deseq_merge$DESeq2.y=="sig up")
sig_agree$DESeq2[idx]=1
idx=which(sigdf_deseq_merge$DESeq2.x=="sig down" & sigdf_deseq_merge$DESeq2.y=="sig down")
sig_agree$DESeq2[idx]=2
sig_agree$edgeR=0
idx=which(sigdf_edgeR_merge$edgeR.x=="sig up" & sigdf_edgeR_merge$edgeR.y=="sig up")
sig_agree$edgeR[idx]=1
idx=which(sigdf_edgeR_merge$edgeR.x=="sig down" & sigdf_edgeR_merge$edgeR.y=="sig down")
sig_agree$edgeR[idx]=2
sig_agree$limmavoom=0
idx=which(sigdf_limmavoom_merge$limmavoom.x=="sig up" & sigdf_limmavoom_merge$limmavoom.y=="sig up")
sig_agree$limmavoom[idx]=1
idx=which(sigdf_limmavoom_merge$limmavoom.x=="sig down" & sigdf_limmavoom_merge$limmavoom.y=="sig down")
sig_agree$limmavoom[idx]=2

sig_df=sig_agree

###########panel B upset plot
library(UpSetR)

list_Input=list('miRglmm up'=rownames(sig_df)[sig_df$miRglmm==1], 'miRglmm down'=rownames(sig_df)[sig_df$miRglmm==2],
                'DESeq2 up'=rownames(sig_df)[sig_df$DESeq2==1], 'DESeq2 down'=rownames(sig_df)[sig_df$DESeq2==2],
                'NB GLM up'=rownames(sig_df)[sig_df$miRglmnb==1], 'NB GLM down'=rownames(sig_df)[sig_df$miRglmnb==2],
                'edgeR up'=rownames(sig_df)[sig_df$edgeR==1], 'edgeR down'=rownames(sig_df)[sig_df$edgeR==2],
                'limma-voom up'=rownames(sig_df)[sig_df$limmavoom==1], 'limma-voom down'=rownames(sig_df)[sig_df$limmavoom==2])
myplot=upset(fromList(list_Input), order.by="freq", sets=rev(names(list_Input)), mainbar.y.label="Number of common significant miRNA",
             sets.x.label="Significant miRNA", text.scale=c(2, 2, 2,2,2,2), keep.order=TRUE)
print(myplot)

upset_ind=myplot[["New_data"]]
x1=unlist(list_Input, use.names=TRUE)
x1=x1[!duplicated(x1)]
idx=which(upset_ind$`miRglmm up`==0 & upset_ind$`miRglmm down`==0 &
            upset_ind$`DESeq2 up`==1 & upset_ind$`DESeq2 down`==0 &
            upset_ind$`edgeR up`==0 & upset_ind$`edgeR down`==0 &
            upset_ind$`limma-voom up`==0 & upset_ind$`limma-voom down`==0 &
            upset_ind$`NB GLM up`==1 & upset_ind$`NB GLM down`==0)
x1[idx]

tiff(file="figures/resub/figure3B.tif", height=9, width=10, units="in", res=320)
print(myplot)
dev.off()


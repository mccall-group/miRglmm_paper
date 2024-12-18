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


miRglmm_kappa=Kappa(table(sigdf_miRglmm_merge$miRglmm.x, sigdf_miRglmm_merge$miRglmm.y))[["Unweighted"]][["value"]]
miRglmnb_kappa=Kappa(table(sigdf_miRglmnb_merge$miRglmnb.x, sigdf_miRglmnb_merge$miRglmnb.y))[["Unweighted"]][["value"]]
deseq_kappa=Kappa(table(sigdf_deseq_merge$DESeq2.x, sigdf_deseq_merge$DESeq2.y))[["Unweighted"]][["value"]]
edgeR_kappa=Kappa(table(sigdf_edgeR_merge$edgeR.x, sigdf_edgeR_merge$edgeR.y))[["Unweighted"]][["value"]]
limmavoom_kappa=Kappa(table(sigdf_limmavoom_merge$limmavoom.x, sigdf_limmavoom_merge$limmavoom.y))[["Unweighted"]][["value"]]



kappa_all=data.frame("kappa"=rbind(miRglmm_kappa, miRglmnb_kappa, deseq_kappa, edgeR_kappa, limmavoom_kappa))
corr_all$kappa=kappa_all[,1]




#create interaction significant groups of all aggregated methods vs miRglmm
sig_groups_v_miRglmm=data.frame("miRglmnb"=as.character(interaction(sig_df$miRglmm, sig_df$miRglmnb)), "DESeq2"=as.character(interaction(sig_df$miRglmm, sig_df$DESeq2)),
                      "edgeR"=as.character(interaction(sig_df$miRglmm, sig_df$edgeR)), "limmavoom"=as.character(interaction(sig_df$miRglmm, sig_df$limmavoom)))
rownames(sig_groups_v_miRglmm)=rownames(sig_df)

sig_groups_v_miRglmm[sig_groups_v_miRglmm=="TRUE.TRUE"]="significant for both"
sig_groups_v_miRglmm[sig_groups_v_miRglmm=="TRUE.FALSE"]="significant for miRglmm only"
sig_groups_v_miRglmm[sig_groups_v_miRglmm=="FALSE.TRUE"]="significant for aggregation method only"
sig_groups_v_miRglmm[sig_groups_v_miRglmm=="FALSE.FALSE"]="significant for neither"

#combine with the betas to make sure miRNAs align with correct beta
comb_df=transform(merge(beta_df, sig_groups_v_miRglmm, by='row.names'), row.names=Row.names, Row.names=NULL)


#calculate ICC and Kappa to include in plotting

icc_out=data.frame(icc=t(sapply(seq(2,5), function(x) icc(comb_df[, c(1,x)], model="oneway", type="agreement", unit="single")$value)))
colnames(icc_out)=colnames(comb_df)[2:5]
kappa_out=data.frame(kappa=t(sapply(seq(2,5), function(x) Kappa(table(sig_df[,1], sig_df[,x]))[["Unweighted"]][["value"]])))
colnames(kappa_out)=colnames(sig_df)[2:5]
#idx=which(rownames(comb_df)=="hsa-miR-27a-3p")
#scatter plots
p1=ggplot(comb_df, aes(x=miRglmnb.x, y=miRglmm, color=miRglmnb.y))+geom_point()+geom_abline()+xlab('NB GLM logFC')+ylab('miRglmm logFC')+scale_colour_discrete(drop=FALSE)+labs(color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out$miRglmnb.x,2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out$miRglmnb,2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  theme(plot.title=element_text(size=15, hjust=0.5), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),
        axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15), legend.position="none")#+
  #annotate("text", x=comb_df$miRglmnb.x[idx], y=comb_df$miRglmm[idx]+0.05, label="o", size=8)


#ggsave("figures/figure4A_legend.tif", plot=last_plot(), device="tiff", scale=6, width=40, height=20, units="mm", dpi=320, bg="white")

p2=ggplot(comb_df, aes(x=DESeq2.x, y=miRglmm, color=DESeq2.y))+geom_point()+geom_abline()+xlab('DESeq2 logFC')+ylab('miRglmm logFC')+scale_colour_discrete(drop=FALSE)+labs(color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out$DESeq2.x,2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out$DESeq2,2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
theme(plot.title=element_text(size=15, hjust=0.5), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),
      axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15), legend.position="none")#+
#annotate("text", x=comb_df$miRglmnb.x[idx], y=comb_df$miRglmm[idx]+0.05, label="o", size=8)


#ggsave("figures/figure4A_2.tif", plot=last_plot(), device="tiff", scale=6, width=40, height=20, units="mm", dpi=320, bg="white")

p3=ggplot(comb_df, aes(x=edgeR.x, y=miRglmm, color=edgeR.y))+geom_point()+geom_abline()+xlab('edgeR logFC')+ylab('miRglmm logFC')+scale_colour_discrete(drop=FALSE)+labs(color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out$edgeR.x,2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out$edgeR,2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  theme(plot.title=element_text(size=15, hjust=0.5), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),
        axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15), legend.position="none")#+
#annotate("text", x=comb_df$miRglmnb.x[idx], y=comb_df$miRglmm[idx]+0.05, label="o", size=8)
  
  #ggsave("figures/figure4A_3.tif", plot=last_plot(), device="tiff", scale=6, width=40, height=20, units="mm", dpi=320, bg="white")

p4=ggplot(comb_df, aes(x=limmavoom.x, y=miRglmm, color=limmavoom.y))+geom_point()+geom_abline()+xlab('limma-voom logFC')+ylab('miRglmm logFC')+scale_colour_discrete(drop=FALSE)+labs(color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out$limmavoom.x,2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out$limmavoom,2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  theme(plot.title=element_text(size=15, hjust=0.5), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),
        axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15), legend.position="none")#+
#annotate("text", x=comb_df$miRglmnb.x[idx], y=comb_df$miRglmm[idx]+0.05, label="o", size=8)
  
p4_legend=ggplot(comb_df, aes(x=limmavoom.x, y=miRglmm, color=limmavoom.y))+geom_point()+geom_abline()+xlab('limma-voom logFC')+ylab('miRglmm logFC')+scale_colour_discrete(drop=FALSE)+labs(color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out$limmavoom.x,2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out$limmavoom,2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  theme(plot.title=element_text(size=15, hjust=0.5), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),
        axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15))#+
#annotate("text", x=comb_df$miRglmnb.x[idx], y=comb_df$miRglmm[idx]+0.05, label="o", size=8)

  ggsave("figures/figure4_legend.tif", plot=last_plot(), device="tiff", scale=6, width=40, height=20, units="mm", dpi=320, bg="white")

###########panel B upset plot
library(UpSetR)

list_Input=list(miRglmm=rownames(sig_df)[sig_df$miRglmm==TRUE], DESeq2=rownames(sig_df)[sig_df$DESeq2==TRUE],`NB GLM`=rownames(sig_df)[sig_df$miRglmnb==TRUE], 
                edgeR=rownames(sig_df)[sig_df$edgeR==TRUE], `limma-voom`=rownames(sig_df)[sig_df$limmavoom==TRUE])
myplot=upset(fromList(list_Input), order.by="freq", mainbar.y.label="Number of common significant miRNA",
             sets.x.label="Significant miRNA", text.scale=c(2, 2, 2,2,2,2))

tiff(file="figures/figure4B.tif", height=6.74, width=8.9, units="in", res=320)
print(myplot)
dev.off()


upset_ind=myplot[["New_data"]]
x1=unlist(list_Input, use.names=TRUE)
x1=x1[!duplicated(x1)]
idx=which(upset_ind$miRglmm==1 & upset_ind$DESeq2==0 & upset_ind$edgeR==0 & upset_ind$`limma-voom`==0 & upset_ind$`NB GLM`==0)
x1[idx]
idx=which(upset_ind$miRglmm==0 & upset_ind$DESeq2==1 & upset_ind$edgeR==1 & upset_ind$`limma-voom`==1 & upset_ind$`NB GLM`==1)
x1[idx]


########## panel C 
load(file='bladder_testes_results/filter_neg1_processed_results.rda')
LRTp=data.frame(results[["LRTp"]])
p5=ggplot(LRTp, aes(x=LRTp))+geom_histogram(color="black", fill="gray", bins=50)+xlab('Likelihood Ratio Test p-value')+ylab('number of miRNA')+
  theme(axis.title.x=element_text(size=15), axis.title.y=element_text(size=15),
        axis.text.x=element_text(size=15), axis.text.y=element_text(size=15))
#ggsave("figures/figure4C.tif", plot=last_plot(), device="tiff", scale=6, width=40, height=20, units="mm", dpi=320, bg="white")


library(ggpubr)
ggarrange(ggarrange(p1, p2, p3, p4, ncol=2, nrow=2), ggarrange(NULL,p5, nrow=2), ncol=2)
ggsave("figures/figure4.tif", plot=last_plot(), device="tiff", width=12.8, height=6.64, units="in", dpi=320, bg="white")

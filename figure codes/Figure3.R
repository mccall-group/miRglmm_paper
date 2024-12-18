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
sig_df=data.frame(padj_df<0.05)
sig_df$contrast=p_df$contrast
colnames(sig_df)=colnames(data.frame(results[["pvals"]]))
rownames(sig_df)=rownames(data.frame(results[["pvals"]]))

beta_df_study89=beta_df
results_study89=results
sig_df_study89=sig_df
padj_df_study89=padj_df


#load miRblood results to compare
load(file='miRblood_exact/results/nofilter_processed_results.rda')

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
sig_df=data.frame(padj_df<0.05)
sig_df$contrast=p_df$contrast
colnames(sig_df)=colnames(data.frame(results[["pvals"]]))
rownames(sig_df)=rownames(data.frame(results[["pvals"]]))

beta_df_miRblood=beta_df
results_miRblood=results
sig_df_miRblood=sig_df
padj_df_miRblood=padj_df

sig_df_miRblood$miRNA=rownames(sig_df_miRblood)
sig_df_study89$miRNA=rownames(sig_df_study89)

#miRglmm study 89 vs miRblood

miRglmm_study89=beta_df_study89[, c("miRglmm" ,"contrast")]
miRglmm_study89$miRNA=rownames(miRglmm_study89)

miRglmm_miRblood=beta_df_miRblood[, c("miRglmm" ,"contrast")]
miRglmm_miRblood$miRNA=rownames(miRglmm_miRblood)

miRglmm_merge=merge(miRglmm_study89, miRglmm_miRblood, by=c("miRNA", "contrast"), all=FALSE)
miRglmm_correlation=miRglmm_merge %>% group_by(contrast) %>% summarise(r=cor(miRglmm.x, miRglmm.y, method="spearman", use="complete.obs"))
miRglmm_correlation$method="miRglmm"
pearson_miRglmm_correlation=miRglmm_merge %>% group_by(contrast) %>% summarise(r=cor(miRglmm.x, miRglmm.y, method="pearson", use="complete.obs"))
pearson_miRglmm_correlation$method="miRglmm"
p1=ggplot(miRglmm_merge, aes(x=miRglmm.x, y=miRglmm.y))+geom_point()+facet_wrap(~contrast, ncol=4, scales="free")+
  theme_minimal()+
  xlab('Juzenas logFC estimate')+ylab('Jehn logFC estimate')+ggtitle('miRglmm')+
  theme(plot.title=element_text(size=12, hjust=0.5))

#NB GLM study 89 vs miRblood

miRglmnb_study89=beta_df_study89[, c("miRglmnb" ,"contrast")]
miRglmnb_study89$miRNA=rownames(miRglmnb_study89)

miRglmnb_miRblood=beta_df_miRblood[, c("miRglmnb" ,"contrast")]
miRglmnb_miRblood$miRNA=rownames(miRglmnb_miRblood)

miRglmnb_merge=merge(miRglmnb_study89, miRglmnb_miRblood, by=c("miRNA", "contrast"))
miRglmnb_correlation=miRglmnb_merge %>% group_by(contrast) %>% summarise(r=cor(miRglmnb.x, miRglmnb.y, method="spearman", use="complete.obs"))
miRglmnb_correlation$method="NB GLM"
pearson_miRglmnb_correlation=miRglmnb_merge %>% group_by(contrast) %>% summarise(r=cor(miRglmnb.x, miRglmnb.y, method="pearson", use="complete.obs"))
pearson_miRglmnb_correlation$method="NB GLM"
p2=ggplot(miRglmnb_merge, aes(x=miRglmnb.x, y=miRglmnb.y))+geom_point()+facet_wrap(~contrast, ncol=4, scales="free")+
  theme_minimal()+
  xlab('Juzenas logFC estimate')+ylab('Jehn logFC estimate')+ggtitle('NB GLM')+
  theme(plot.title=element_text(size=12, hjust=0.5))

#deseq2 study 89 vs miRblood

deseq_study89=beta_df_study89[, c("DESeq2" ,"contrast")]
deseq_study89$miRNA=rownames(deseq_study89)

deseq_miRblood=beta_df_miRblood[, c("DESeq2" ,"contrast")]
deseq_miRblood$miRNA=rownames(deseq_miRblood)

deseq_merge=merge(deseq_study89, deseq_miRblood, by=c("miRNA", "contrast"))
deseq_correlation=deseq_merge %>% group_by(contrast) %>% summarise(r=cor(DESeq2.x, DESeq2.y, method="spearman", use="complete.obs"))
deseq_correlation$method="DESeq2"
pearson_deseq_correlation=deseq_merge %>% group_by(contrast) %>% summarise(r=cor(DESeq2.x, DESeq2.y, method="pearson", use="complete.obs"))
pearson_deseq_correlation$method="DESeq2"
p3=ggplot(deseq_merge, aes(x=DESeq2.x, y=DESeq2.y))+geom_point()+facet_wrap(~contrast, ncol=4, scales="free")+
  theme_minimal()+
  xlab('Juzenas logFC estimate')+ylab('Jehn logFC estimate')+ggtitle('DESEq2')+
  theme(plot.title=element_text(size=12, hjust=0.5))

#edgeR study 89 vs miRblood

edgeR_study89=beta_df_study89[, c("edgeR" ,"contrast")]
edgeR_study89$miRNA=rownames(edgeR_study89)

edgeR_miRblood=beta_df_miRblood[, c("edgeR" ,"contrast")]
edgeR_miRblood$miRNA=rownames(edgeR_miRblood)

edgeR_merge=merge(edgeR_study89, edgeR_miRblood, by=c("miRNA", "contrast"))
edgeR_correlation=edgeR_merge %>% group_by(contrast) %>% summarise(r=cor(edgeR.x, edgeR.y, method="spearman", use="complete.obs"))
edgeR_correlation$method="edgeR"
pearson_edgeR_correlation=edgeR_merge %>% group_by(contrast) %>% summarise(r=cor(edgeR.x, edgeR.y, method="pearson", use="complete.obs"))
pearson_edgeR_correlation$method="edgeR"
p4=ggplot(edgeR_merge, aes(x=edgeR.x, y=edgeR.y))+geom_point()+facet_wrap(~contrast, ncol=4, scales="free")+
  theme_minimal()+
  xlab('Juzenas logFC estimate')+ylab('Jehn logFC estimate')+ggtitle('edgeR')+
  theme(plot.title=element_text(size=12, hjust=0.5))

#limma-voom study 89 vs miRblood

limmavoom_study89=beta_df_study89[, c("limmavoom" ,"contrast")]
limmavoom_study89$miRNA=rownames(limmavoom_study89)

limmavoom_miRblood=beta_df_miRblood[, c("limmavoom" ,"contrast")]
limmavoom_miRblood$miRNA=rownames(limmavoom_miRblood)

limmavoom_merge=merge(limmavoom_study89, limmavoom_miRblood, by=c("miRNA", "contrast"))
limmavoom_correlation=limmavoom_merge %>% group_by(contrast) %>% summarise(r=cor(limmavoom.x, limmavoom.y, method="spearman", use="complete.obs"))
limmavoom_correlation$method="limma-voom"
pearson_limmavoom_correlation=limmavoom_merge %>% group_by(contrast) %>% summarise(r=cor(limmavoom.x, limmavoom.y, method="pearson", use="complete.obs"))
pearson_limmavoom_correlation$method="limma-voom"
p5=ggplot(limmavoom_merge, aes(x=limmavoom.x, y=limmavoom.y))+geom_point()+facet_wrap(~contrast, ncol=4, scales="free")+
  theme_minimal()+
  xlab('Juzenas logFC estimate')+ylab('Jehn logFC estimate')+ggtitle('limma-voom')+
  theme(plot.title=element_text(size=12, hjust=0.5))

library(ggpubr)
ggarrange(p1, p2, p3, p4, p5, ncol=1)

corr_all=rbind(miRglmm_correlation, miRglmnb_correlation, deseq_correlation, edgeR_correlation, limmavoom_correlation)
corr_mat=data.frame(spread(corr_all, key=method, value=r))
corr_mat=corr_mat[, c("contrast", "miRglmm", "NB.GLM", "DESeq2", "edgeR", "limma.voom")]
colnames(corr_mat)=c("contrast", "miRglmm", "NB GLM", "DESeq2", "edgeR", "limma-voom")

pearson_corr_all=rbind(pearson_miRglmm_correlation, pearson_miRglmnb_correlation, pearson_deseq_correlation, pearson_edgeR_correlation, pearson_limmavoom_correlation)
pearson_corr_mat=data.frame(spread(pearson_corr_all, key=method, value=r))
pearson_corr_mat=pearson_corr_mat[, c("contrast", "miRglmm", "NB.GLM", "DESeq2", "edgeR", "limma.voom")]
colnames(pearson_corr_mat)=c("contrast", "miRglmm", "NB GLM", "DESeq2", "edgeR", "limma-voom")


sigdf_miRglmm_merge=merge(sig_df_study89[, c("miRglmm" ,"contrast", "miRNA")], sig_df_miRblood[, c("miRglmm" ,"contrast", "miRNA")], by=c("miRNA", "contrast"))
sigdf_miRglmnb_merge=merge(sig_df_study89[, c("miRglmnb" ,"contrast", "miRNA")], sig_df_miRblood[, c("miRglmnb" ,"contrast", "miRNA")], by=c("miRNA", "contrast"))
sigdf_deseq_merge=merge(sig_df_study89[, c("DESeq2" ,"contrast", "miRNA")], sig_df_miRblood[, c("DESeq2" ,"contrast", "miRNA")], by=c("miRNA", "contrast"))
sigdf_edgeR_merge=merge(sig_df_study89[, c("edgeR" ,"contrast", "miRNA")], sig_df_miRblood[, c("edgeR" ,"contrast", "miRNA")], by=c("miRNA", "contrast"))
sigdf_limmavoom_merge=merge(sig_df_study89[, c("limmavoom" ,"contrast", "miRNA")], sig_df_miRblood[, c("limmavoom" ,"contrast", "miRNA")], by=c("miRNA", "contrast"))


miRglmm_kappa=sigdf_miRglmm_merge %>% group_by(contrast) %>% summarise(kappa=Kappa(table(miRglmm.x, miRglmm.y))[["Unweighted"]][["value"]])
miRglmnb_kappa=sigdf_miRglmnb_merge %>% group_by(contrast) %>% summarise(kappa=Kappa(table(miRglmnb.x, miRglmnb.y))[["Unweighted"]][["value"]])
deseq_kappa=sigdf_deseq_merge %>% group_by(contrast) %>% summarise(kappa=Kappa(table(DESeq2.x, DESeq2.y))[["Unweighted"]][["value"]])
edgeR_kappa=sigdf_edgeR_merge %>% group_by(contrast) %>% summarise(kappa=Kappa(table(edgeR.x, edgeR.y))[["Unweighted"]][["value"]])
limmavoom_kappa=sigdf_limmavoom_merge %>% group_by(contrast) %>% summarise(kappa=Kappa(table(limmavoom.x, limmavoom.y))[["Unweighted"]][["value"]])

miRglmm_kappa$method="miRglmm"
miRglmnb_kappa$method="NB GLM"
deseq_kappa$method="DESeq2"
edgeR_kappa$method="edgeR"
limmavoom_kappa$method="limma-voom"

kappa_all=rbind(miRglmm_kappa, miRglmnb_kappa, deseq_kappa, edgeR_kappa, limmavoom_kappa)
kappa_mat=data.frame(spread(kappa_all, key=method, value=kappa))
kappa_mat=kappa_mat[, c("contrast", "miRglmm", "NB.GLM", "DESeq2", "edgeR", "limma.voom")]
colnames(kappa_mat)=c("contrast", "miRglmm", "NB GLM", "DESeq2", "edgeR", "limma-voom")


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


#calculate ICC and Kappa to include in plotting

icc_out=data.frame(icc=t(sapply(seq(2,5), function(x) icc(comb_df[, c(1,x)], model="oneway", type="agreement", unit="single")$value)))
colnames(icc_out)=colnames(comb_df)[2:5]
kappa_out=data.frame(kappa=t(sapply(seq(2,5), function(x) Kappa(table(sig_df[,1], sig_df[,x]))[["Unweighted"]][["value"]])))
colnames(kappa_out)=colnames(sig_df)[2:5]



#idx=which(rownames(comb_df)=="hsa-miR-27a-3p")
#scatter plots
p1=ggplot(comb_df, aes(x=miRglmnb.x, y=miRglmm, color=miRglmnb.y, shape=contrast.x))+geom_point(size=3)+geom_abline()+xlab('NB GLM logFC')+ylab('miRglmm logFC')+scale_colour_discrete(drop=FALSE)+labs(color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out$miRglmnb.x,2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out$miRglmnb,2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  theme(plot.title=element_text(size=15, hjust=0.5), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),
        axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15), legend.position="none")#+
  #annotate("text", x=comb_df$miRglmnb.x[idx], y=comb_df$miRglmm[idx]+0.05, label="o", size=8)


#ggsave("figures/figure4A_legend.tif", plot=last_plot(), device="tiff", scale=6, width=40, height=20, units="mm", dpi=320, bg="white")

p2=ggplot(comb_df, aes(x=DESeq2.x, y=miRglmm, color=DESeq2.y, shape=contrast.x))+geom_point(size=3)+geom_abline()+xlab('DESeq2 logFC')+ylab('miRglmm logFC')+scale_colour_discrete(drop=FALSE)+labs(color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out$DESeq2.x,2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out$DESeq2,2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
theme(plot.title=element_text(size=15, hjust=0.5), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),
      axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15), legend.position="none")#+
#annotate("text", x=comb_df$miRglmnb.x[idx], y=comb_df$miRglmm[idx]+0.05, label="o", size=8)


#ggsave("figures/figure4A_2.tif", plot=last_plot(), device="tiff", scale=6, width=40, height=20, units="mm", dpi=320, bg="white")

p3=ggplot(comb_df, aes(x=edgeR.x, y=miRglmm, color=edgeR.y, shape=contrast.x))+geom_point(size=3)+geom_abline()+xlab('edgeR logFC')+ylab('miRglmm logFC')+scale_colour_discrete(drop=FALSE)+labs(color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out$edgeR.x,2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out$edgeR,2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  theme(plot.title=element_text(size=15, hjust=0.5), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),
        axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15), legend.position="none")#+
#annotate("text", x=comb_df$miRglmnb.x[idx], y=comb_df$miRglmm[idx]+0.05, label="o", size=8)
  
  #ggsave("figures/figure4A_3.tif", plot=last_plot(), device="tiff", scale=6, width=40, height=20, units="mm", dpi=320, bg="white")

p4=ggplot(comb_df, aes(x=limmavoom.x, y=miRglmm, color=limmavoom.y, shape=contrast.x))+geom_point(size=3)+geom_abline()+xlab('limma-voom logFC')+ylab('miRglmm logFC')+scale_colour_discrete(drop=FALSE)+labs(color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out$limmavoom.x,2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out$limmavoom,2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  theme(plot.title=element_text(size=15, hjust=0.5), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),
        axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15), legend.position="none")#+
#annotate("text", x=comb_df$miRglmnb.x[idx], y=comb_df$miRglmm[idx]+0.05, label="o", size=8)
  
p4_legend=ggplot(comb_df, aes(x=limmavoom.x, y=miRglmm, color=limmavoom.y, shape=contrast.x))+geom_point(size=3)+geom_abline()+xlab('limma-voom logFC')+ylab('miRglmm logFC')+scale_colour_discrete(drop=FALSE)+labs(color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out$limmavoom.x,2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out$limmavoom,2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  theme(plot.title=element_text(size=15, hjust=0.5), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),
        axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15))#+
#annotate("text", x=comb_df$miRglmnb.x[idx], y=comb_df$miRglmm[idx]+0.05, label="o", size=8)

  ggsave("figures/figure3_legend.tif", plot=last_plot(), device="tiff", scale=6, width=40, height=20, units="mm", dpi=320, bg="white")

###########panel B upset plot
library(UpSetR)

#find edgeR significance from fit object
#load(file='bladder_testes_results/filter_neg1_results.rda')
#sig_df_edgeR=data.frame('edgeR'=fits[["edgeR"]][["table"]][["PValue"]]<0.05)
#rownames(sig_df_edgeR)=rownames(fits[["edgeR"]][["table"]])
#sig_df=transform(merge(sig_df, sig_df_edgeR, by='row.names'), row.names=Row.names, Row.names=NULL)




list_Input=list(miRglmm=rownames(sig_df)[sig_df$miRglmm==TRUE], DESeq2=rownames(sig_df)[sig_df$DESeq2==TRUE],`NB GLM`=rownames(sig_df)[sig_df$miRglmnb==TRUE], 
                edgeR=rownames(sig_df)[sig_df$edgeR==TRUE], `limma-voom`=rownames(sig_df)[sig_df$limmavoom==TRUE])

myplot=upset(fromList(list_Input), order.by="freq", mainbar.y.label="Number of common significant miRNA",
             sets.x.label="Significant miRNA", text.scale=c(2, 2, 2,2,2,2))

tiff(file="figures/figure3B.tif", height=6.74, width=8.9, units="in", res=320)
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
LRTp=data.frame(results[["LRTp"]])
p5=ggplot(LRTp, aes(x=LRTp))+geom_histogram(color="black", fill="gray", bins=50)+xlab('Likelihood Ratio Test p-value')+ylab('number of miRNA')+
  theme(axis.title.x=element_text(size=15), axis.title.y=element_text(size=15),
        axis.text.x=element_text(size=15), axis.text.y=element_text(size=15))
#ggsave("figures/figure4C.tif", plot=last_plot(), device="tiff", scale=6, width=40, height=20, units="mm", dpi=320, bg="white")


library(ggpubr)
ggarrange(ggarrange(p1, p2, p3, p4, ncol=2, nrow=2), ggarrange(NULL,p5, nrow=2), ncol=2)
ggsave("figures/figure3.tif", plot=last_plot(), device="tiff", width=12.8, height=6.64, units="in", dpi=320, bg="white")

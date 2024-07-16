library(irr)
library(vcd)
library(ggplot2)

load(file='study 89/results/filter_neg1_processed_results.rda')
contrast_in="Monocyte vs B lymphocyte"

beta_df=data.frame(results[["beta_hat"]])
beta_df=beta_df[which(beta_df$contrast==contrast_in),]
p_df=data.frame(results[["pvals"]])
p_df=p_df[which(p_df$contrast==contrast_in),]
#adjust pvalue for number of miRNA analyzed (dont correct for all contrasts performed)
padj_df=matrix(, nrow=dim(p_df)[1], ncol=dim(p_df)[2]-1)
for (ind in seq(1,dim(p_df)[2]-1)){
  padj_df[,ind]=p.adjust(p_df[,ind], method="BH")
 }
sig_df=data.frame(padj_df<0.05)
sig_df$contrast=p_df$contrast
colnames(sig_df)=colnames(data.frame(p_df))
rownames(sig_df)=rownames(data.frame(p_df))




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
p1=ggplot(comb_df, aes(x=miRglmnb.x, y=miRglmm, color=miRglmnb.y))+geom_point(size=3)+geom_abline()+xlab('NB GLM logFC')+ylab('miRglmm logFC')+scale_colour_discrete(drop=FALSE)+labs(color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out$miRglmnb.x,2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out$miRglmnb,2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  theme(plot.title=element_text(size=15, hjust=0.5), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),
        axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15), legend.position="none")#+
  #annotate("text", x=comb_df$miRglmnb.x[idx], y=comb_df$miRglmm[idx]+0.05, label="o", size=8)



p2=ggplot(comb_df, aes(x=DESeq2.x, y=miRglmm, color=DESeq2.y))+geom_point(size=3)+geom_abline()+xlab('DESeq2 logFC')+ylab('miRglmm logFC')+scale_colour_discrete(drop=FALSE)+labs(color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out$DESeq2.x,2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out$DESeq2,2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
theme(plot.title=element_text(size=15, hjust=0.5), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),
      axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15), legend.position="none")#+
#annotate("text", x=comb_df$miRglmnb.x[idx], y=comb_df$miRglmm[idx]+0.05, label="o", size=8)



p3=ggplot(comb_df, aes(x=edgeR.x, y=miRglmm, color=edgeR.y))+geom_point(size=3)+geom_abline()+xlab('edgeR logFC')+ylab('miRglmm logFC')+scale_colour_discrete(drop=FALSE)+labs(color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out$edgeR.x,2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out$edgeR,2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  theme(plot.title=element_text(size=15, hjust=0.5), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),
        axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15), legend.position="none")#+
#annotate("text", x=comb_df$miRglmnb.x[idx], y=comb_df$miRglmm[idx]+0.05, label="o", size=8)
  

p4=ggplot(comb_df, aes(x=limmavoom.x, y=miRglmm, color=limmavoom.y))+geom_point(size=3)+geom_abline()+xlab('limma-voom logFC')+ylab('miRglmm logFC')+scale_colour_discrete(drop=FALSE)+labs(color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out$limmavoom.x,2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out$limmavoom,2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  theme(plot.title=element_text(size=15, hjust=0.5), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),
        axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15), legend.position="none")#+
#annotate("text", x=comb_df$miRglmnb.x[idx], y=comb_df$miRglmm[idx]+0.05, label="o", size=8)
  
p4_legend=ggplot(comb_df, aes(x=limmavoom.x, y=miRglmm, color=limmavoom.y))+geom_point(size=3)+geom_abline()+xlab('limma-voom logFC')+ylab('miRglmm logFC')+scale_colour_discrete(drop=FALSE)+labs(color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out$limmavoom.x,2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out$limmavoom,2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  theme(plot.title=element_text(size=15, hjust=0.5), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),
        axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15))#+
#annotate("text", x=comb_df$miRglmnb.x[idx], y=comb_df$miRglmm[idx]+0.05, label="o", size=8)

  ggsave("figures/supfigure6_legend.tif", plot=last_plot(), device="tiff", scale=6, width=40, height=20, units="mm", dpi=320, bg="white")

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

tiff(file="figures/supfigure6B.tif", height=6.74, width=8.9, units="in", res=320)
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
ggsave("figures/supfigure6.tif", plot=last_plot(), device="tiff", width=12.8, height=6.64, units="in", dpi=320, bg="white")

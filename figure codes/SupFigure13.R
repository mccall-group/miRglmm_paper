library(irr)
library(vcd)

load(file='bladder_testes_results/filter_neg1_processed_results.rda')

beta_df=data.frame(results[["beta_hat"]])
p_df=data.frame(results[["pvals"]])
padj_df=matrix(, nrow=dim(p_df)[1], ncol=dim(p_df)[2])
for (ind in seq(1,dim(p_df)[2])){
  padj_df[,ind]=p.adjust(p_df[,ind], method="BH")
}
sig_df=data.frame(padj_df<0.05)
colnames(sig_df)=colnames(data.frame(results[["pvals"]]))
rownames(sig_df)=rownames(data.frame(results[["pvals"]]))


#create interaction significant groups of all aggregated methods vs miRglmm
sig_groups_v_miRglmm=data.frame("int1"=as.character(interaction(sig_df$DESeq2, sig_df$miRglmnb)), "int2"=as.character(interaction(sig_df$edgeR, sig_df$miRglmnb)),
                      "int3"=as.character(interaction(sig_df$limmavoom, sig_df$miRglmnb)), "int4"=as.character(interaction(sig_df$edgeR, sig_df$DESeq2)),
                      "int5"=as.character(interaction(sig_df$limmavoom, sig_df$DESeq2)), "int6"=as.character(interaction(sig_df$limmavoom, sig_df$edgeR)))
rownames(sig_groups_v_miRglmm)=rownames(sig_df)

sig_groups_v_miRglmm[sig_groups_v_miRglmm=="TRUE.TRUE"]="significant for both"
sig_groups_v_miRglmm[sig_groups_v_miRglmm=="TRUE.FALSE"]="significant for y-axis method only"
sig_groups_v_miRglmm[sig_groups_v_miRglmm=="FALSE.TRUE"]="significant for x-axis method only"
sig_groups_v_miRglmm[sig_groups_v_miRglmm=="FALSE.FALSE"]="significant for neither"

#combine with the betas to make sure miRNAs align with correct beta
comb_df=transform(merge(beta_df, sig_groups_v_miRglmm, by='row.names'), row.names=Row.names, Row.names=NULL)

comb_df$int1=factor(comb_df$int1, levels=c("significant for x-axis method only", "significant for both",
                                           "significant for y-axis method only", "significant for neither"))
comb_df$int2=factor(comb_df$int2, levels=c("significant for x-axis method only", "significant for both",
                                           "significant for y-axis method only", "significant for neither"))

comb_df$int3=factor(comb_df$int3, levels=c("significant for x-axis method only", "significant for both",
                                           "significant for y-axis method only", "significant for neither"))

comb_df$int4=factor(comb_df$int4, levels=c("significant for x-axis method only", "significant for both",
                                           "significant for y-axis method only", "significant for neither"))

comb_df$int5=factor(comb_df$int5, levels=c("significant for x-axis method only", "significant for both",
                                           "significant for y-axis method only", "significant for neither"))

comb_df$int6=factor(comb_df$int6, levels=c("significant for x-axis method only", "significant for both",
                                           "significant for y-axis method only", "significant for neither"))

#calculate ICC and Kappa to include in plotting

icc_out=data.frame(icc=t(sapply(seq(3,5), function(x) icc(comb_df[, c(2,x)], model="oneway", type="agreement", unit="single")$value)))
colnames(icc_out)=colnames(comb_df)[3:5]
icc_out2=data.frame(icc=t(sapply(seq(4,5), function(x) icc(comb_df[, c(3,x)], model="oneway", type="agreement", unit="single")$value)))
colnames(icc_out2)=colnames(comb_df)[4:5]
icc_out3=data.frame(icc=t(sapply(5, function(x) icc(comb_df[, c(4,x)], model="oneway", type="agreement", unit="single")$value)))
colnames(icc_out3)=colnames(comb_df)[5]
icc_out=cbind(icc_out, icc_out2, icc_out3)

kappa_out=data.frame(kappa=t(sapply(seq(3,5), function(x) Kappa(table(sig_df[,2], sig_df[,x]))[["Unweighted"]][["value"]])))
colnames(kappa_out)=colnames(sig_df)[3:5]
kappa_out2=data.frame(kappa=t(sapply(seq(4,5), function(x) Kappa(table(sig_df[,3], sig_df[,x]))[["Unweighted"]][["value"]])))
colnames(kappa_out2)=colnames(sig_df)[4:5]
kappa_out3=data.frame(kappa=t(sapply(5, function(x) Kappa(table(sig_df[,4], sig_df[,x]))[["Unweighted"]][["value"]])))
colnames(kappa_out3)=colnames(sig_df)[5]
kappa_out=cbind(kappa_out, kappa_out2, kappa_out3)
#idx=which(rownames(comb_df)=="hsa-miR-27a-3p")
#scatter plots
p1=ggplot(comb_df, aes(x=miRglmnb, y=DESeq2, color=int1))+geom_point()+geom_abline()+xlab('NB GLM logFC')+ylab('DESeq2 logFC')+scale_colour_discrete(drop=FALSE)+labs(title="NB GLM vs DESeq2", color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out[1],2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out[1],2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  theme(plot.title=element_text(size=15, hjust=0.5), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),
        axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.position="none", 
        legend.title=element_text(size=15), legend.text=element_text(size=15))#+
  #annotate("text", x=comb_df$miRglmnb.x[idx], y=comb_df$miRglmm[idx]+0.05, label="o", size=8)


#ggsave("figures/supfigure4_1.tif", plot=last_plot(), device="tiff", scale=6, width=40, height=20, units="mm", dpi=320, bg="white")

p2=ggplot(comb_df, aes(x=miRglmnb, y=edgeR, color=int2))+geom_point()+geom_abline()+xlab('NB GLM logFC')+ylab('edgeR logFC')+scale_colour_discrete(drop=FALSE)+labs(title="NB GLM vs edgeR", color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out[2],2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out[2],2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  theme(plot.title=element_text(size=15, hjust=0.5), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),
        axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.position="none", 
        legend.title=element_text(size=15), legend.text=element_text(size=15))#+
#annotate("text", x=comb_df$miRglmnb.x[idx], y=comb_df$miRglmm[idx]+0.05, label="o", size=8)


#ggsave("figures/supfigure4_2.tif", plot=last_plot(), device="tiff", scale=6, width=40, height=20, units="mm", dpi=320, bg="white")


p3=ggplot(comb_df, aes(x=miRglmnb, y=limmavoom, color=int3))+geom_point()+geom_abline()+xlab('NB GLM logFC')+ylab('limma-voom logFC')+scale_colour_discrete(drop=FALSE)+labs(title="NB GLM vs limma-voom", color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out[3],2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out[3],2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  theme(plot.title=element_text(size=15, hjust=0.5), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),
        axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.position="none", 
        legend.title=element_text(size=15), legend.text=element_text(size=15))#+
#annotate("text", x=comb_df$miRglmnb.x[idx], y=comb_df$miRglmm[idx]+0.05, label="o", size=8)


#ggsave("figures/supfigure4_3.tif", plot=last_plot(), device="tiff", scale=6, width=40, height=20, units="mm", dpi=320, bg="white")


p4=ggplot(comb_df, aes(x=DESeq2, y=edgeR, color=int4))+geom_point()+geom_abline()+xlab('DESeq2 logFC')+ylab('edgeR logFC')+scale_colour_discrete(drop=FALSE)+labs(title="DESeq2 vs edgeR", color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out[4],2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out[4],2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  theme(plot.title=element_text(size=15, hjust=0.5), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),
        axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.position="none", 
        legend.title=element_text(size=15), legend.text=element_text(size=15))#+
#annotate("text", x=comb_df$miRglmnb.x[idx], y=comb_df$miRglmm[idx]+0.05, label="o", size=8)


#ggsave("figures/supfigure4_4.tif", plot=last_plot(), device="tiff", scale=6, width=40, height=20, units="mm", dpi=320, bg="white")

p5=ggplot(comb_df, aes(x=DESeq2, y=limmavoom, color=int5))+geom_point()+geom_abline()+xlab('DESeq2 logFC')+ylab('limma-voom logFC')+scale_colour_discrete(drop=FALSE)+labs(title="DESeq2 vs limma-voom", color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out[5],2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out[5],2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  theme(plot.title=element_text(size=15, hjust=0.5), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),
        axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.position="none", 
        legend.title=element_text(size=15), legend.text=element_text(size=15))#+
#annotate("text", x=comb_df$miRglmnb.x[idx], y=comb_df$miRglmm[idx]+0.05, label="o", size=8)


#ggsave("figures/supfigure4_5.tif", plot=last_plot(), device="tiff", scale=6, width=40, height=20, units="mm", dpi=320, bg="white")

p6=ggplot(comb_df, aes(x=edgeR, y=limmavoom, color=int6))+geom_point()+geom_abline()+xlab('edgeR logFC')+ylab('limma-voom logFC')+scale_colour_discrete(drop=FALSE)+labs(title="edgeR vs limma-voom", color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out[6],2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out[6],2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  theme(plot.title=element_text(size=15, hjust=0.5), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),
        axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.position="none", 
        legend.title=element_text(size=15), legend.text=element_text(size=15))#+
#annotate("text", x=comb_df$miRglmnb.x[idx], y=comb_df$miRglmm[idx]+0.05, label="o", size=8)

ggplot(comb_df, aes(x=edgeR, y=limmavoom, color=int6))+geom_point()+geom_abline()+xlab('edgeR logFC')+ylab('limma-voom logFC')+scale_colour_discrete(drop=FALSE)+labs(title="edgeR vs limma-voom", color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out[6],2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out[6],2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)+
  theme(plot.title=element_text(size=15, hjust=0.5), axis.text.x=element_text(size=15), axis.title.x=element_text(size=15),
        axis.text.y=element_text(size=15), axis.title.y=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15))#+
#annotate("text", x=comb_df$miRglmnb.x[idx], y=comb_df$miRglmm[idx]+0.05, label="o", size=8)
ggsave("figures/supfigure13_legend.tif", plot=last_plot(), device="tiff", width=6.98, height=4.84, units="in", dpi=320, bg="white")

library(ggpubr)
ggarrange(p1, p2, p3, p4, p5, p6, nrow=3, ncol=2)

ggsave("figures/supfigure13.tif", plot=last_plot(), device="tiff", width=12, height=14, units="in", dpi=320, bg="white")


### summarize simulation results by truth via boxplots
library(reshape2)
library(ggplot2)


############# MSE by truth boxplots
load(file="sim results/sims_N100_m2_s1_rtruncnorm13_combinedresults.rda")
uniq_miRNA=seq(1,length(results[["MSE_sim_by_truth"]]))

out_miRglmm=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["MSE_sim_by_truth"]][[row]][["miRglmm"]]))))
colnames(out_miRglmm)=results[["MSE_sim_by_truth"]][[1]][["true_beta"]]
out_miRglmm$method="miRglmm"

out_miRglmnb=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["MSE_sim_by_truth"]][[row]][["miRglmnb"]]))))
colnames(out_miRglmnb)=results[["MSE_sim_by_truth"]][[1]][["true_beta"]]
out_miRglmnb$method="NB GLM"

out_DESeq2=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["MSE_sim_by_truth"]][[row]][["DESeq2"]]))))
colnames(out_DESeq2)=results[["MSE_sim_by_truth"]][[1]][["true_beta"]]
out_DESeq2$method="DESeq2"

out_edgeR=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["MSE_sim_by_truth"]][[row]][["edgeR"]]))))
colnames(out_edgeR)=results[["MSE_sim_by_truth"]][[1]][["true_beta"]]
out_edgeR$method="edgeR"

out_limvoom=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["MSE_sim_by_truth"]][[row]][["limmavoom"]]))))
colnames(out_limvoom)=results[["MSE_sim_by_truth"]][[1]][["true_beta"]]
out_limvoom$method="limma-voom"

out_all=rbind(out_miRglmm, out_miRglmnb, out_DESeq2, out_edgeR, out_limvoom)
data_long = melt(out_all, 
                 id.vars = c("method"),
                 variable.name = "truth", 
                 value.name = "MSE")

data_long$`True LogFC`=as.factor(round(exp(as.numeric(as.character(data_long$truth))), digits=2))
data_long$method=factor(data_long$method, levels=c("miRglmm", "empty", "NB GLM", "DESeq2", "edgeR", "limma-voom"))
data_long$MSE=data_long$MSE*1000
p1=ggplot(data_long, aes(x=`True LogFC`, y=MSE, fill=method))+geom_boxplot()+ylab(expression(paste("MSE (10"^"-3", ")")))+
  xlab('True Fold Change')+
  scale_fill_discrete(drop=FALSE, labels=c("miRglmm", "NB GLM", "DESeq2", "edgeR", "limma-voom"), 
                      breaks=c("miRglmm", "NB GLM", "DESeq2", "edgeR", "limma-voom"))+
                                        theme(axis.title=element_text(size=15),
                                         axis.text=element_text(size=15),
                                         legend.title=element_text(size=15),
                                         legend.text=element_text(size=15),
                                         legend.position="bottom")
print(p1)
#ggsave("figures/figure15A.tif", plot=last_plot(), device="tiff", width=300, height=150, units="mm", dpi=320, bg="white")


############# coverage probability by truth boxplots

uniq_miRNA=seq(1,length(results[["coverage_prob_sim_by_truth"]]))

out_miRglmm=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["coverage_prob_sim_by_truth"]][[row]][["miRglmm"]]))))
colnames(out_miRglmm)=results[["coverage_prob_sim_by_truth"]][[1]][["true_beta"]]
out_miRglmm$method="miRglmm"

out_miRglmnb=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["coverage_prob_sim_by_truth"]][[row]][["miRglmnb"]]))))
colnames(out_miRglmnb)=results[["coverage_prob_sim_by_truth"]][[1]][["true_beta"]]
out_miRglmnb$method="NB GLM"

out_DESeq2=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["coverage_prob_sim_by_truth"]][[row]][["DESeq2"]]))))
colnames(out_DESeq2)=results[["coverage_prob_sim_by_truth"]][[1]][["true_beta"]]
out_DESeq2$method="DESeq2"

out_limvoom=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["coverage_prob_sim_by_truth"]][[row]][["limmavoom"]]))))
colnames(out_limvoom)=results[["coverage_prob_sim_by_truth"]][[1]][["true_beta"]]
out_limvoom$method="limma-voom"

out_all=rbind(out_miRglmm, out_miRglmnb, out_DESeq2, out_limvoom)
data_long = melt(out_all, 
                 id.vars = c("method"),
                 variable.name = "truth", 
                 value.name = "coverage probability")

data_long$`True LogFC`=as.factor(round(exp(as.numeric(as.character(data_long$truth))), digits=3))
data_long$method=factor(data_long$method, levels=c("miRglmm", "empty", "NB GLM", "DESeq2", "edgeR", "limma-voom"))
p2=ggplot(data_long, aes(x=`True LogFC`, y=`coverage probability`, fill=method))+geom_boxplot()+
  ylab('coverage proportion')+xlab('True Fold Change')+
  scale_fill_discrete(drop=FALSE, labels=c("miRglmm", "NB GLM", "DESeq2", "limma-voom"), 
                      breaks=c("miRglmm", "NB GLM", "DESeq2", "limma-voom"))+
  theme(axis.title=element_text(size=15), axis.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        legend.position="bottom")
print(p2)
#ggsave("figures/figure15B.tif", plot=last_plot(), device="tiff", width=300, height=150, units="mm", dpi=320, bg="white")


############# proportion of significant random effects by truth
uniq_miRNA=seq(1,length(results[["prop_sig_by_truth"]]))

out_all=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["prop_sig_by_truth"]][[row]][["sig"]]))))
colnames(out_all)=results[["prop_sig_by_truth"]][[1]][["true_beta"]]
data_long = melt(out_all, 
                 variable.name = "truth", 
                 value.name = "proportion significant")
data_long$`True LogFC`=as.factor(round(exp(as.numeric(as.character(data_long$truth))), digits=3))

p3=ggplot(data_long, aes(x=`True LogFC`, y=`proportion significant`))+geom_boxplot()+
  ylab(paste('proportion of miRNA with', '\n' , 'significant random slope effect'))+xlab('True Fold Change')+
  theme(axis.title=element_text(size=15), axis.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        legend.position="bottom")
print(p3)
library(ggpubr)
#ggarrange(p1, p2, p3, nrow=2, ncol=2)
ggarrange(p1, p2, nrow=1, ncol=2)
ggsave("figures/figure2.tif", plot=last_plot(), device="tiff", width=17, height=11, units="in", dpi=320, bg="white")

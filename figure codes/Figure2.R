### summarize simulation results by truth via boxplots
library(reshape2)
library(ggplot2)
library(dplyr)

############# MSE by truth boxplots
load(file="sim results/sims_N100_m2_s1_rtruncnorm13_combinedresults.rda")
uniq_miRNA=seq(1,length(results[["MSE_sim_by_truth"]]))

out_miRglmm=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["MSE_sim_by_truth"]][[row]][["miRglmm"]]))))
colnames(out_miRglmm)=results[["MSE_sim_by_truth"]][[1]][["true_beta"]]
out_miRglmm$method="miRglmm"

out_miRglmm_pois=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["MSE_sim_by_truth"]][[row]][["miRglmm_poiss"]]))))
colnames(out_miRglmm_pois)=results[["MSE_sim_by_truth"]][[1]][["true_beta"]]
out_miRglmm_pois$method="miRglmm-Poisson"

out_miRglmnb=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["MSE_sim_by_truth"]][[row]][["miRglmnb"]]))))
colnames(out_miRglmnb)=results[["MSE_sim_by_truth"]][[1]][["true_beta"]]
out_miRglmnb$method="NB GLM"

# out_miRglmpois=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["MSE_sim_by_truth"]][[row]][["miRglmpois"]]))))
# colnames(out_miRglmpois)=results[["MSE_sim_by_truth"]][[1]][["true_beta"]]
# out_miRglmpois$method="Poisson GLM"

out_DESeq2=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["MSE_sim_by_truth"]][[row]][["DESeq2"]]))))
colnames(out_DESeq2)=results[["MSE_sim_by_truth"]][[1]][["true_beta"]]
out_DESeq2$method="DESeq2"

out_edgeR=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["MSE_sim_by_truth"]][[row]][["edgeR"]]))))
colnames(out_edgeR)=results[["MSE_sim_by_truth"]][[1]][["true_beta"]]
out_edgeR$method="edgeR"

out_limvoom=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["MSE_sim_by_truth"]][[row]][["limmavoom"]]))))
colnames(out_limvoom)=results[["MSE_sim_by_truth"]][[1]][["true_beta"]]
out_limvoom$method="limma-voom"

out_all=rbind(out_miRglmm, out_miRglmm_pois, out_miRglmnb, 
              #out_miRglmpois, 
              out_DESeq2, out_edgeR, out_limvoom)
data_long = melt(out_all, 
                 id.vars = c("method"),
                 variable.name = "truth", 
                 value.name = "MSE")

data_long$`True LogFC`=as.factor(round(exp(as.numeric(as.character(data_long$truth))), digits=2))
data_long$method=factor(data_long$method, levels=c("miRglmm", "empty", "miRglmm-Poisson", "NB GLM", "DESeq2", "edgeR", "limma-voom", "empty2"))
data_long$MSE=data_long$MSE*1000
p1=ggplot(data_long, aes(x=`True LogFC`, y=MSE, fill=method))+geom_boxplot()+ylab(expression(paste("MSE (10"^"-3", ")")))+
  xlab('True Fold Change')+theme_minimal()+
  scale_fill_discrete(drop=FALSE, labels=c("miRglmm", "miRglmm-Poisson", "NB GLM", "DESeq2", "edgeR", "limma-voom"), 
                      breaks=c("miRglmm", "miRglmm-Poisson", "NB GLM", "DESeq2", "edgeR", "limma-voom"))+
                                        theme(axis.title=element_text(size=15),
                                         axis.text=element_text(size=15),
                                         legend.title=element_text(size=15),
                                         legend.text=element_text(size=15),
                                         legend.position="none")
print(p1)
#ggsave("figures/figure15A.tif", plot=last_plot(), device="tiff", width=300, height=150, units="mm", dpi=320, bg="white")


############# coverage probability by truth boxplots

uniq_miRNA=seq(1,length(results[["coverage_prob_sim_by_truth"]]))

out_miRglmm=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["coverage_prob_sim_by_truth"]][[row]][["miRglmm"]]))))
colnames(out_miRglmm)=results[["coverage_prob_sim_by_truth"]][[1]][["true_beta"]]
out_miRglmm$method="miRglmm"

out_miRglmm_pois=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["coverage_prob_sim_by_truth"]][[row]][["miRglmm_poiss"]]))))
colnames(out_miRglmm_pois)=results[["coverage_prob_sim_by_truth"]][[1]][["true_beta"]]
out_miRglmm_pois$method="miRglmm-Poisson"

out_miRglmnb=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["coverage_prob_sim_by_truth"]][[row]][["miRglmnb"]]))))
colnames(out_miRglmnb)=results[["coverage_prob_sim_by_truth"]][[1]][["true_beta"]]
out_miRglmnb$method="NB GLM"

# out_miRglmpois=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["coverage_prob_sim_by_truth"]][[row]][["miRglmpois"]]))))
# colnames(out_miRglmpois)=results[["coverage_prob_sim_by_truth"]][[1]][["true_beta"]]
# out_miRglmpois$method="Poisson GLM"

out_DESeq2=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["coverage_prob_sim_by_truth"]][[row]][["DESeq2"]]))))
colnames(out_DESeq2)=results[["coverage_prob_sim_by_truth"]][[1]][["true_beta"]]
out_DESeq2$method="DESeq2"

out_limvoom=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["coverage_prob_sim_by_truth"]][[row]][["limmavoom"]]))))
colnames(out_limvoom)=results[["coverage_prob_sim_by_truth"]][[1]][["true_beta"]]
out_limvoom$method="limma-voom"

out_all=rbind(out_miRglmm, out_miRglmm_pois, out_miRglmnb, 
             # out_miRglmpois,
              out_DESeq2, out_limvoom)
data_long = melt(out_all, 
                 id.vars = c("method"),
                 variable.name = "truth", 
                 value.name = "coverage probability")

data_long$`True LogFC`=as.factor(round(exp(as.numeric(as.character(data_long$truth))), digits=3))
data_long$method=factor(data_long$method, levels=c("miRglmm", "empty", "miRglmm-Poisson", "NB GLM", "DESeq2", "edgeR", "limma-voom", "empty2"))
p2=ggplot(data_long, aes(x=`True LogFC`, y=`coverage probability`, fill=method))+geom_boxplot()+
  ylab('coverage proportion')+xlab('True Fold Change')+theme_minimal()+
  scale_fill_discrete(drop=FALSE, labels=c("miRglmm", "miRglmm-Poisson", "NB GLM", "DESeq2", "limma-voom"), 
                      breaks=c("miRglmm", "miRglmm-Poisson", "NB GLM", "DESeq2", "limma-voom"))+
  theme(axis.title=element_text(size=15), axis.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        legend.position="none")
print(p2)
#ggsave("figures/figure15B.tif", plot=last_plot(), device="tiff", width=300, height=150, units="mm", dpi=320, bg="white")

### variance of estimates by truth
beta_var_by_truth=lapply(results[["beta_hat"]], function(x) x %>% group_by(true_beta) %>% summarise_all(funs(var)))
uniq_miRNA=seq(1,length(beta_var_by_truth))

out_miRglmm=data.frame(rbind(t(sapply(uniq_miRNA, function(row) beta_var_by_truth[[row]][["miRglmm"]]))))
colnames(out_miRglmm)=beta_var_by_truth[[1]][["true_beta"]]
out_miRglmm$method="miRglmm"

out_miRglmm_pois=data.frame(rbind(t(sapply(uniq_miRNA, function(row) beta_var_by_truth[[row]][["miRglmm_poiss"]]))))
colnames(out_miRglmm_pois)=beta_var_by_truth[[1]][["true_beta"]]
out_miRglmm_pois$method="miRglmm-Poisson"

out_miRglmnb=data.frame(rbind(t(sapply(uniq_miRNA, function(row) beta_var_by_truth[[row]][["miRglmnb"]]))))
colnames(out_miRglmnb)=beta_var_by_truth[[1]][["true_beta"]]
out_miRglmnb$method="NB GLM"

# out_miRglmpois=data.frame(rbind(t(sapply(uniq_miRNA, function(row) beta_var_by_truth[[row]][["miRglmpois"]]))))
# colnames(out_miRglmpois)=beta_var_by_truth[[1]][["true_beta"]]
# out_miRglmpois$method="Poisson GLM"

out_DESeq2=data.frame(rbind(t(sapply(uniq_miRNA, function(row) beta_var_by_truth[[row]][["DESeq2"]]))))
colnames(out_DESeq2)=beta_var_by_truth[[1]][["true_beta"]]
out_DESeq2$method="DESeq2"

out_limvoom=data.frame(rbind(t(sapply(uniq_miRNA, function(row) beta_var_by_truth[[row]][["limmavoom"]]))))
colnames(out_limvoom)=beta_var_by_truth[[1]][["true_beta"]]
out_limvoom$method="limma-voom"

out_all=rbind(out_miRglmm, out_miRglmm_pois, out_miRglmnb, 
              #out_miRglmpois,
              out_DESeq2, out_limvoom)
data_long = melt(out_all, 
                 id.vars = c("method"),
                 variable.name = "truth", 
                 value.name = "variance")
data_long$variance=data_long$variance*1000

data_long$`True LogFC`=as.factor(round(exp(as.numeric(as.character(data_long$truth))), digits=3))
data_long$method=factor(data_long$method, levels=c("miRglmm", "empty", "miRglmm-Poisson", "NB GLM", "DESeq2", "edgeR", "limma-voom", "empty2"))
p3=ggplot(data_long, aes(x=`True LogFC`, y=`variance`, fill=method))+geom_boxplot()+
  ylab(expression(paste("variance (10"^"-3", ")")))+xlab('True Fold Change')+theme_minimal()+
  scale_fill_discrete(drop=FALSE, labels=c("miRglmm", "miRglmm-Poisson", "NB GLM", "DESeq2", "limma-voom"), 
                      breaks=c("miRglmm", "miRglmm-Poisson", "NB GLM", "DESeq2", "limma-voom"))+
  theme(axis.title=element_text(size=15), axis.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        legend.position="none")
print(p3)

# ### variance of estimates by truth
# CIwidth_by_truth=lapply(results[["CI_width"]], function(x) cbind(x, "true_beta"=results[["beta_hat"]][[1]]$true_beta) %>% group_by(true_beta) %>% summarise_all(funs(mean)))
# uniq_miRNA=seq(1,length(CIwidth_by_truth))
# 
# out_miRglmm=data.frame(rbind(t(sapply(uniq_miRNA, function(row) CIwidth_by_truth[[row]][["miRglmm"]]))))
# colnames(out_miRglmm)=CIwidth_by_truth[[1]][["true_beta"]]
# out_miRglmm$method="miRglmm"
# 
# out_miRglmm_pois=data.frame(rbind(t(sapply(uniq_miRNA, function(row) CIwidth_by_truth[[row]][["miRglmm_poiss"]]))))
# colnames(out_miRglmm_pois)=CIwidth_by_truth[[1]][["true_beta"]]
# out_miRglmm_pois$method="miRglmm-Poisson"
# 
# out_miRglmnb=data.frame(rbind(t(sapply(uniq_miRNA, function(row) CIwidth_by_truth[[row]][["miRglmnb"]]))))
# colnames(out_miRglmnb)=CIwidth_by_truth[[1]][["true_beta"]]
# out_miRglmnb$method="NB GLM"
# 
# out_miRglmpois=data.frame(rbind(t(sapply(uniq_miRNA, function(row) CIwidth_by_truth[[row]][["miRglmpois"]]))))
# colnames(out_miRglmpois)=CIwidth_by_truth[[1]][["true_beta"]]
# out_miRglmpois$method="Poisson GLM"
# 
# out_DESeq2=data.frame(rbind(t(sapply(uniq_miRNA, function(row) CIwidth_by_truth[[row]][["DESeq2"]]))))
# colnames(out_DESeq2)=CIwidth_by_truth[[1]][["true_beta"]]
# out_DESeq2$method="DESeq2"
# 
# out_limvoom=data.frame(rbind(t(sapply(uniq_miRNA, function(row) CIwidth_by_truth[[row]][["limmavoom"]]))))
# colnames(out_limvoom)=CIwidth_by_truth[[1]][["true_beta"]]
# out_limvoom$method="limma-voom"
# 
# out_all=rbind(out_miRglmm, out_miRglmm_pois, out_miRglmnb, out_miRglmpois,out_DESeq2, out_limvoom)
# data_long = melt(out_all, 
#                  id.vars = c("method"),
#                  variable.name = "truth", 
#                  value.name = "CI width")
# 
# data_long$`True LogFC`=as.factor(round(exp(as.numeric(as.character(data_long$truth))), digits=3))
# data_long$method=factor(data_long$method, levels=c("miRglmm", "empty", "miRglmm-Poisson", "NB GLM", "Poisson GLM", "DESeq2", "edgeR", "limma-voom", "empty2"))
# p4=ggplot(data_long, aes(x=`True LogFC`, y=`CI width`, fill=method))+geom_boxplot()+
#   ylab('CI width')+xlab('True Fold Change')+
#   scale_fill_discrete(drop=FALSE, labels=c("miRglmm", "miRglmm-Poisson", "NB GLM", "Poisson GLM", "DESeq2", "limma-voom"), 
#                       breaks=c("miRglmm", "miRglmm-Poisson", "NB GLM", "Poisson GLM", "DESeq2", "limma-voom"))+
#   theme(axis.title=element_text(size=15), axis.text=element_text(size=15),
#         legend.title=element_text(size=15),
#         legend.text=element_text(size=15),
#         legend.position="bottom")
# print(p4)


#####proportion of significant miRNA
uniq_miRNA=seq(1,length(results[["sig_by_truth"]]))
results[["FDR"]]=lapply(results[["pvals"]], function(x) data.frame(matrix(unlist(lapply(x, function(y) p.adjust(y, method="fdr"))), ncol=9)))



TPR_pos=lapply(results[["FDR"]], function(x) t(data.frame(TPR=colSums(x[which(x[,9]>0),]<0.05)/length(which(x[,9]>0)))))
TPR_neg=lapply(results[["FDR"]], function(x) t(data.frame(TPR=colSums(x[which(x[,9]<0),]<0.05)/length(which(x[,9]<0)))))
TNR_sim=lapply(results[["FDR"]], function(x) t(data.frame(TNR=colSums(x[which(x[,9]==0),]>0.05)/length(which(x[,9]==0)))))

TPR_pos_mat=as.data.frame(matrix(unlist(TPR_pos), ncol=9, byrow=TRUE))
colnames(TPR_pos_mat)=colnames(results[["pvals"]][[1]])
TPR_pos_mat=TPR_pos_mat[, 1:8]

TPR_neg_mat=as.data.frame(matrix(unlist(TPR_neg), ncol=9, byrow=TRUE))
colnames(TPR_neg_mat)=colnames(results[["pvals"]][[1]])
TPR_neg_mat=TPR_neg_mat[, 1:8]

TNR_mat=as.data.frame(matrix(unlist(TNR_sim), ncol=9, byrow=TRUE))
colnames(TNR_mat)=colnames(results[["pvals"]][[1]])
TNR_mat=TNR_mat[, 1:8]

out_miRglmm=data.frame(cbind(TPR_neg_mat$miRglmm, TNR_mat$miRglmm, TPR_pos_mat$miRglmm))
colnames(out_miRglmm)=results[["sig_by_truth"]][[1]][["true_beta"]]
out_miRglmm$method="miRglmm"

out_miRglmm_pois=data.frame(cbind(TPR_neg_mat$miRglmm_poiss, TNR_mat$miRglmm_poiss, TPR_pos_mat$miRglmm_poiss))
colnames(out_miRglmm_pois)=results[["sig_by_truth"]][[1]][["true_beta"]]
out_miRglmm_pois$method="miRglmm-Poisson"

out_miRglmnb=data.frame(cbind(TPR_neg_mat$miRglmnb, TNR_mat$miRglmnb, TPR_pos_mat$miRglmnb))
colnames(out_miRglmnb)=results[["sig_by_truth"]][[1]][["true_beta"]]
out_miRglmnb$method="NB GLM"

# out_miRglmpois=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["sig_by_truth"]][[row]][["miRglmpois"]]))))
# colnames(out_miRglmpois)=results[["sig_by_truth"]][[1]][["true_beta"]]
# out_miRglmpois$method="Poisson GLM"

out_DESeq2=data.frame(cbind(TPR_neg_mat$DESeq2, TNR_mat$DESeq2, TPR_pos_mat$DESeq2))
colnames(out_DESeq2)=results[["sig_by_truth"]][[1]][["true_beta"]]
out_DESeq2$method="DESeq2"

out_limvoom=data.frame(cbind(TPR_neg_mat$limmavoom, TNR_mat$limmavoom, TPR_pos_mat$limmavoom))
colnames(out_limvoom)=results[["sig_by_truth"]][[1]][["true_beta"]]
out_limvoom$method="limma-voom"

out_wilcox=data.frame(cbind(TPR_neg_mat$wilcoxp, TNR_mat$wilcoxp, TPR_pos_mat$wilcoxp))
colnames(out_wilcox)=results[["sig_by_truth"]][[1]][["true_beta"]]
out_wilcox$method="Wilcoxon"

out_all=rbind(out_miRglmm, out_miRglmm_pois, out_miRglmnb, 
              #out_miRglmpois,
              out_DESeq2, out_limvoom, out_wilcox)
data_long = melt(out_all, 
                 id.vars = c("method"),
                 variable.name = "truth", 
                 value.name = "proportion of significant miRNA")

data_long$`True LogFC`=as.factor(round(exp(as.numeric(as.character(data_long$truth))), digits=3))
data_long$method=factor(data_long$method, levels=c("miRglmm", "empty", "miRglmm-Poisson", "NB GLM", "DESeq2", "edgeR", "limma-voom", "Wilcoxon"))
p4=ggplot(data_long, aes(x=`True LogFC`, y=`proportion of significant miRNA`, fill=method))+geom_boxplot()+
  ylab('proportion of significant miRNA')+xlab('True Fold Change')+theme_minimal()+
  scale_fill_discrete(drop=FALSE, labels=c("miRglmm", "miRglmm-Poisson", "NB GLM", "DESeq2", "limma-voom", "Wilcoxon"), 
                      breaks=c("miRglmm", "miRglmm-Poisson", "NB GLM", "DESeq2", "limma-voom", "Wilcoxon"))+
  theme(axis.title=element_text(size=15), axis.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        legend.position="none")
print(p4)
p4_legend=ggplot(data_long, aes(x=`True LogFC`, y=`proportion of significant miRNA`, fill=method))+geom_boxplot()+
  ylab('proportion of significant miRNA')+xlab('True Fold Change')+theme_minimal()+
  scale_fill_discrete(drop=FALSE, labels=c("miRglmm", "miRglmm-Poisson", "NB GLM", "DESeq2", "limma-voom", "Wilcoxon"), 
                      breaks=c("miRglmm", "miRglmm-Poisson", "NB GLM", "DESeq2", "limma-voom", "Wilcoxon"))+
  theme(axis.title=element_text(size=15), axis.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        legend.position="bottom")
# ############# proportion of significant random effects by truth
# uniq_miRNA=seq(1,length(results[["prop_sig_by_truth"]]))
# 
# out_all=data.frame(rbind(t(sapply(uniq_miRNA, function(row) results[["prop_sig_by_truth"]][[row]][["sig"]]))))
# colnames(out_all)=results[["prop_sig_by_truth"]][[1]][["true_beta"]]
# data_long = melt(out_all, 
#                  variable.name = "truth", 
#                  value.name = "proportion significant")
# data_long$`True LogFC`=as.factor(round(exp(as.numeric(as.character(data_long$truth))), digits=3))
# 
# p3=ggplot(data_long, aes(x=`True LogFC`, y=`proportion significant`))+geom_boxplot()+
#   ylab(paste('proportion of miRNA with', '\n' , 'significant random slope effect'))+xlab('True Fold Change')+
#   theme(axis.title=element_text(size=15), axis.text=element_text(size=15),
#         legend.title=element_text(size=15),
#         legend.text=element_text(size=15),
#         legend.position="bottom")
# print(p3)
library(ggpubr)
#ggarrange(p1, p2, p3, nrow=2, ncol=2)
ggarrange(p1, p2,p3,p4, nrow=2, ncol=2)
ggsave("figures/figure2.tif", plot=last_plot(), device="tiff", width=17, height=11, units="in", dpi=320, bg="white")
print(p4_legend)
ggsave("figures/figure2_legend.tif", plot=last_plot(), device="tiff", width=17, height=11, units="in", dpi=320, bg="white")

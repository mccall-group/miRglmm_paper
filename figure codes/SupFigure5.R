library(ggplot2)

########boxplots
#load truth data
ratio_pool=read.csv(file="ERCC files/FINAL_Ratiometric_SynthA_and_SynthB-1_fixlast3.csv", sep="\t")
ratio_pool$A=as.numeric(gsub("x","",ratio_pool$X.SAMPLE.A.))
ratio_pool$B=as.numeric(gsub("x","",ratio_pool$X.SAMPLE.B.))
ratio_pool$ratio=ratio_pool$B/ratio_pool$A
ratio_pool$true_logFC=log(ratio_pool$ratio)
out=ratio_pool["true_logFC"]
rownames(out)=ratio_pool$ratio.seqID



#load processed results
load('ERCC/allfilters_processed_results.rda')


beta_hat=results[["beta_hat"]]
beta_hat=transform(merge(beta_hat, out, by='row.names'), row.names=Row.names, Row.names=NULL)
beta_hat=beta_hat[,c("filter..1", "no.filter","miRglmnb","DESeq2","edgeR","limmavoom", "true_logFC")]
colnames(beta_hat)=c("miRglmm filter -1","miRglmm no filter","NB GLM","DESeq2","edgeR","limma-voom", "True LogFC")

diff=beta_hat[,!(colnames(beta_hat)=="True LogFC")]-beta_hat[,"True LogFC"]
diff_long=stack(diff)
colnames(diff_long)=c("diff", "method")
diff_long$`True LogFC`=as.factor(round(as.numeric(as.character(rep(exp(beta_hat[,"True LogFC"]), times=ncol(diff)))), digits=2))

p1=ggplot(diff_long, aes(x=`True LogFC`, y=diff, fill=method))+geom_boxplot()+
  ylab('Pool LogFC Estimate-Truth')+xlab('True Pool Fold Change')+theme(axis.title=element_text(size=15),
                                                                              axis.text=element_text(size=15),
                                                                              legend.title=element_text(size=15),
                                                                              legend.text=element_text(size=15))
p1_nolegend=ggplot(diff_long, aes(x=`True LogFC`, y=diff, fill=method))+geom_boxplot()+
  ylab('Pool LogFC Estimate-Truth')+xlab('True Pool Fold Change')+theme(axis.title=element_text(size=15),
                                                                         axis.text=element_text(size=15),
                                                                         legend.position="none")

#ggsave("figures/figure6A.tif", plot=last_plot(), device="tiff", width=8.4, height=4.88, units="in", dpi=320, bg="white")

######## MSE by filter
library(tidyverse)
library(reshape2)
library(lme4)
library(SummarizedExperiment)
library(DESeq2)
library(edgeR)

#load truth data
ratio_pool=read.csv(file="ERCC files/FINAL_Ratiometric_SynthA_and_SynthB-1_fixlast3.csv", sep="\t")
ratio_pool$A=as.numeric(gsub("x","",ratio_pool$X.SAMPLE.A.))
ratio_pool$B=as.numeric(gsub("x","",ratio_pool$X.SAMPLE.B.))
ratio_pool$ratio=ratio_pool$B/ratio_pool$A
ratio_pool$true_logFC=log(ratio_pool$ratio)
out=ratio_pool["true_logFC"]
rownames(out)=ratio_pool$ratio.seqID



#load processed results
load('ERCC/filter_agg_betahats.rda')

filters=c("no filter", "-1", "-0.5", "0", "0.5", "1", "1.5", "2")


#combine filtered agg with filtered seq models to get table of betas for each model
for (ind in seq(1, length(filters))){
  filter_in=filters[ind]
  beta_hat_in=beta_hat[[filter_in]]
  beta_hat_in=transform(merge(beta_hat_in, out, by='row.names'), row.names=Row.names, Row.names=NULL)
  diff_squared=(beta_hat_in[,!(colnames(beta_hat_in)=="true_logFC")]-beta_hat_in[,"true_logFC"])^2
  
  if (ind==1){
    MSE=data.frame(colMeans(diff_squared, na.rm=TRUE))
  } else {
    MSE=cbind(MSE, data.frame(colMeans(diff_squared, na.rm=TRUE)))
  }
  
}
idx=which(rownames(MSE)=="miRglmnb")
rownames(MSE)[idx]="NB GLM"
idx=which(rownames(MSE)=="limmavoom")
rownames(MSE)[idx]="limma-voom"
colnames(MSE)=filters
MSE_long=stack(MSE)
colnames(MSE_long)=c("MSE", "filter")
MSE_long$method=rep(rownames(MSE), times=ncol(MSE))
MSE_long$method=factor(MSE_long$method, levels=c("miRglmm", "miRglmm2", "NB GLM", "DESeq2", "edgeR", "limma-voom"))
MSE_long$MSE=MSE_long$MSE*1000
p2=ggplot(MSE_long, aes(x=filter, y=MSE, group=method, color=method))+
  scale_color_discrete(drop=FALSE,labels=c("miRglmm", "NB GLM", "DESeq2", "edgeR", "limma-voom"), 
                                  breaks=c("miRglmm", "NB GLM", "DESeq2", "edgeR", "limma-voom"))+
  geom_line(size=2)+xlab('filter (log(median CPM))')+theme(axis.title=element_text(size=15),
                    axis.text=element_text(size=15),
                    legend.title=element_text(size=15),
                    legend.text=element_text(size=15),
                    legend.position="bottom")+ylab(expression(paste("MSE (10"^"-3", ")")))

p2_nolegend=ggplot(MSE_long, aes(x=filter, y=MSE, group=method, color=method))+
  scale_color_discrete(drop=FALSE, labels=c("miRglmm", "NB GLM", "DESeq2", "limma-voom"), 
                                   breaks=c("miRglmm", "NB GLM", "DESeq2", "limma-voom"))+
  geom_line(size=2)+xlab('filter (log(median CPM))')+theme(axis.title=element_text(size=15),
                          axis.text=element_text(size=15),
                          legend.position="none")+ylab(expression(paste("MSE (10"^"-3", ")")))

#ggsave("figures/figure6B.tif", plot=last_plot(), device="tiff", width=4.2, height=4.88, units="in", dpi=320, bg="white")

######## histogram
library(ggplot2)
load(file='ERCC/allfilters_processed_results.rda')
LRTp=data.frame(LRTp=results[["LRTp"]][["filter..1"]])
p3=ggplot(LRTp, aes(x=LRTp))+geom_histogram(color="black", fill="gray", bins=50)+
  xlab('Likelihood Ratio Test p-value')+ylab('number of miRNA')+xlim(c(0,NA))+theme(axis.title=element_text(size=15),
                                                             axis.text=element_text(size=15),
                                                             legend.title=element_text(size=15),
                                                             legend.text=element_text(size=15))



#ggsave("figures/figure6C.tif", plot=last_plot(), device="tiff", width=4.2, height=4.88, units="in", dpi=320, bg="white")

library(ggpubr)
#ggarrange(p1, ggarrange(p2,p3,ncol=2), nrow=2)
plot(p1)
ggsave("figures/figure5A_legend.tif", plot=last_plot(), device="tiff", width=11.83, height=3.25, units="in", dpi=320, bg="white")

plot(p2)
ggsave("figures/figure5B_legend.tif", plot=last_plot(), device="tiff", width=11.83, height=3.25, units="in", dpi=320, bg="white")


ggarrange(p1_nolegend, ggarrange(p2_nolegend,p3,ncol=2), nrow=2)
ggsave("figures/SupFigure5.tif", plot=last_plot(), device="tiff", width=11.83, height=6.53, units="in", dpi=320, bg="white")

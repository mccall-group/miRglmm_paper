library(ggplot2)
library(lme4)
library(scales)


## read in results
load('bladder_testes_results/filter_neg1_processed_results.rda')
beta_hat=results[["beta_hat"]]

############################# significant for miRglmm only
miRNA_plot="hsa-miR-100-5p"#"hsa-miR-100-5p"#

#extract testes FC
miRNA_in=as.data.frame(exp(unlist(beta_hat[which(rownames(beta_hat)==miRNA_plot),])))
miRNA_in$method=rownames(miRNA_in)
colnames(miRNA_in)=c("estimate", "method")
miRNA_in$Pool=rep("Testes", dim(miRNA_in)[1])
#set Bladder FC to 1
miRNA_in_bladder=miRNA_in
miRNA_in_bladder$Pool=rep("Bladder", dim(miRNA_in)[1])
miRNA_in_bladder$estimate=rep(1, dim(miRNA_in)[1])


comb_df=rbind(miRNA_in, miRNA_in_bladder)
comb_df$method=factor(comb_df$method, levels=c("miRglmm", "miRglmm2", "miRglmnb", "DESeq2", "edgeR", "limmavoom"))

ggplot(comb_df, aes(x=Pool, y=estimate, group=method, color=method))+geom_line(size=2)+
  scale_color_discrete(drop=FALSE,labels=c("miRglmm", "NB GLM", "DESeq2", "edgeR", "limma-voom"), 
                                breaks=c("miRglmm", "miRglmnb", "DESeq2", "edgeR", "limmavoom"))+
  ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.05,0))+xlab('')+
 theme(legend.position="bottom", legend.direction="horizontal", 
        plot.title=element_text(size=20, hjust=0.5), axis.title=element_text(size=20), axis.text=element_text(size=15), legend.text=element_text(size=15))

ggsave("figures/supfigure14_legend.tif", plot=last_plot(), device="tiff", width=12.8, height=4.3, units="in", dpi=320, bg="white")
  

p1=ggplot(comb_df, aes(x=Pool, y=estimate, group=method, color=method))+geom_line(size=2)+
  scale_color_discrete(drop=FALSE,labels=c("miRglmm", "NB GLM", "DESeq2", "edgeR", "limma-voom"), 
                       breaks=c("miRglmm", "miRglmnb", "DESeq2", "edgeR", "limmavoom"))+
  scale_y_continuous(breaks=c(0.2, 0.4, 0.6, 0.8, 1.0), labels=c(0.2, 0.4, 0.6, 0.8, 1.0))+
  ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.05,0))+xlab('')+
  ggtitle("hsa-miR-100-5p")+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))


############################# significant for miRglmm only
miRNA_plot="hsa-miR-143-5p"#"hsa-miR-100-5p"#

#extract testes FC
miRNA_in=as.data.frame(exp(unlist(beta_hat[which(rownames(beta_hat)==miRNA_plot),])))
miRNA_in$method=rownames(miRNA_in)
colnames(miRNA_in)=c("estimate", "method")
miRNA_in$Pool=rep("Testes", dim(miRNA_in)[1])
#set Bladder FC to 1
miRNA_in_bladder=miRNA_in
miRNA_in_bladder$Pool=rep("Bladder", dim(miRNA_in)[1])
miRNA_in_bladder$estimate=rep(1, dim(miRNA_in)[1])


comb_df=rbind(miRNA_in, miRNA_in_bladder)
comb_df$method=factor(comb_df$method, levels=c("miRglmm", "miRglmm2", "miRglmnb", "DESeq2", "edgeR", "limmavoom"))

p2=ggplot(comb_df, aes(x=Pool, y=estimate, group=method, color=method))+geom_line(size=2)+
  scale_color_discrete(drop=FALSE,labels=c("miRglmm", "NB GLM", "DESeq2", "edgeR", "limma-voom"), 
                       breaks=c("miRglmm", "miRglmnb", "DESeq2", "edgeR", "limmavoom"))+
  ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.05,0))+xlab('')+
  ggtitle("hsa-miR-143-5p")+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))

############################# significant for miRglmm only
miRNA_plot="hsa-miR-222-3p"#"hsa-miR-100-5p"#

#extract testes FC
miRNA_in=as.data.frame(exp(unlist(beta_hat[which(rownames(beta_hat)==miRNA_plot),])))
miRNA_in$method=rownames(miRNA_in)
colnames(miRNA_in)=c("estimate", "method")
miRNA_in$Pool=rep("Testes", dim(miRNA_in)[1])
#set Bladder FC to 1
miRNA_in_bladder=miRNA_in
miRNA_in_bladder$Pool=rep("Bladder", dim(miRNA_in)[1])
miRNA_in_bladder$estimate=rep(1, dim(miRNA_in)[1])


comb_df=rbind(miRNA_in, miRNA_in_bladder)
comb_df$method=factor(comb_df$method, levels=c("miRglmm", "miRglmm2", "miRglmnb", "DESeq2", "edgeR", "limmavoom"))

data_sub=comb_df[comb_df$method=="miRglmm",]

p3=ggplot(comb_df, aes(x=Pool, y=estimate, group=method, color=method))+geom_line(size=2)+
  scale_color_discrete(drop=FALSE,labels=c("miRglmm", "NB GLM", "DESeq2", "edgeR", "limma-voom"), 
                       breaks=c("miRglmm", "miRglmnb", "DESeq2", "edgeR", "limmavoom"))+
  ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.05,0))+xlab('')+
  scale_y_continuous(breaks=c(1, 2, 3, 4), labels=label_number(accuracy=0.1))+
  ggtitle("hsa-miR-222-3p")+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(), axis.title=element_text(size=20), 
        axis.text=element_text(size=20))+
        geom_line(data_sub, mapping=aes(x=Pool, y=estimate, group=method, color=method), size=2)


  

############################# miRglmm only method not sig

miRNA_plot="hsa-miR-25-3p"#"hsa-miR-100-5p"#

#extract testes FC
miRNA_in=as.data.frame(exp(unlist(beta_hat[which(rownames(beta_hat)==miRNA_plot),])))
miRNA_in$method=rownames(miRNA_in)
colnames(miRNA_in)=c("estimate", "method")
miRNA_in$Pool=rep("Testes", dim(miRNA_in)[1])
#set Bladder FC to 1
miRNA_in_bladder=miRNA_in
miRNA_in_bladder$Pool=rep("Bladder", dim(miRNA_in)[1])
miRNA_in_bladder$estimate=rep(1, dim(miRNA_in)[1])


comb_df=rbind(miRNA_in, miRNA_in_bladder)
comb_df$method=factor(comb_df$method, levels=c("miRglmm", "miRglmm2", "miRglmnb", "DESeq2", "edgeR", "limmavoom"))



p4=ggplot(comb_df, aes(x=Pool, y=estimate, group=method, color=method))+geom_line(size=2)+
  scale_color_discrete(drop=FALSE,labels=c("miRglmm", "NB GLM", "DESeq2", "edgeR", "limma-voom"), 
                       breaks=c("miRglmm", "miRglmnb", "DESeq2", "edgeR", "limmavoom"))+
  scale_y_continuous(breaks=c(2,4,6,8), labels=label_number(accuracy=0.1))+
  ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.05,0))+xlab('')+
  ggtitle("hsa-miR-25-3p")+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))


############################# miRglmm only method not sig
miRNA_plot="hsa-miR-423-5p"#

#extract testes FC
miRNA_in=as.data.frame(exp(unlist(beta_hat[which(rownames(beta_hat)==miRNA_plot),])))
miRNA_in$method=rownames(miRNA_in)
colnames(miRNA_in)=c("estimate", "method")
miRNA_in$Pool=rep("Testes", dim(miRNA_in)[1])
#set Bladder FC to 1
miRNA_in_bladder=miRNA_in
miRNA_in_bladder$Pool=rep("Bladder", dim(miRNA_in)[1])
miRNA_in_bladder$estimate=rep(1, dim(miRNA_in)[1])


comb_df=rbind(miRNA_in, miRNA_in_bladder)
comb_df$method=factor(comb_df$method, levels=c("miRglmm", "miRglmm2", "miRglmnb", "DESeq2", "edgeR", "limmavoom"))


p5=ggplot(comb_df, aes(x=Pool, y=estimate, group=method, color=method))+geom_line(size=2)+
  scale_color_discrete(drop=FALSE,labels=c("miRglmm", "NB GLM", "DESeq2", "edgeR", "limma-voom"), 
                       breaks=c("miRglmm", "miRglmnb", "DESeq2", "edgeR", "limmavoom"))+
  scale_y_continuous(breaks=c(2,4,6), labels=label_number(accuracy=0.1))+
  ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.05,0))+xlab('')+
  ggtitle("hsa-miR-423-5p")+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5),
        axis.title.x=element_blank(), axis.title=element_text(size=20), axis.text=element_text(size=20))

############################# miRglmm only method not sig
miRNA_plot="hsa-miR-664a-5p.SNPC"#

#extract testes FC
miRNA_in=as.data.frame(exp(unlist(beta_hat[which(rownames(beta_hat)==miRNA_plot),])))
miRNA_in$method=rownames(miRNA_in)
colnames(miRNA_in)=c("estimate", "method")
miRNA_in$Pool=rep("Testes", dim(miRNA_in)[1])
#set Bladder FC to 1
miRNA_in_bladder=miRNA_in
miRNA_in_bladder$Pool=rep("Bladder", dim(miRNA_in)[1])
miRNA_in_bladder$estimate=rep(1, dim(miRNA_in)[1])


comb_df=rbind(miRNA_in, miRNA_in_bladder)
comb_df$method=factor(comb_df$method, levels=c("miRglmm", "miRglmm2", "miRglmnb", "DESeq2", "edgeR", "limmavoom"))



p6=ggplot(comb_df, aes(x=Pool, y=estimate, group=method, color=method))+geom_line(size=2)+
  scale_color_discrete(drop=FALSE,labels=c("miRglmm", "NB GLM", "DESeq2", "edgeR", "limma-voom"), 
                       breaks=c("miRglmm", "miRglmnb", "DESeq2", "edgeR", "limmavoom"))+
  scale_y_continuous(breaks=c(5, 10, 15, 20))+
  ylab(paste('Expression relative to', '\n', 'average bladder expression'))+scale_x_discrete(expand=c(0.05,0))+xlab('')+
  ggtitle("hsa-miR-664a-5p")+
  theme(legend.position="none",
        plot.title=element_text(size=20, hjust=0.5), 
        axis.title.x=element_blank(),
        axis.title=element_text(size=20), axis.text=element_text(size=20))



library(ggpubr)
ggarrange(p1, p2, p3, p4, p5, p6, nrow=2, ncol=3)

ggsave("figures/supfigure14.tif", plot=last_plot(), device="tiff", width=26, height=12, units="in", dpi=320, bg="white")

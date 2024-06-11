library(vegan)
library(ggplot2)
library(SummarizedExperiment)
########## NMDS plot

load("ERCC/ERCC_filtered.rda")


Bray_curtis_NMDS=metaMDS(as.matrix(t(assay(panel_B_filter))), k=2)
out=Bray_curtis_NMDS$points
col_data=as.data.frame(colData(panel_B_filter))
col_data$NMDS1=out[,1]
col_data$NMDS2=out[,2]


ggplot(col_data, aes(x=NMDS1, y=NMDS2, color=Lab))+geom_point(aes(shape=Pool, size=4))+
  scale_size(guide='none')+guides(color=guide_legend(override.aes=list(size=3)))+
  guides(shape=guide_legend(override.aes=list(size=3)))+theme(legend.position="bottom", legend.text=element_text(size=15), axis.title=element_text(size=15))

#ggsave("figures/supfigure3_A.tif", plot=last_plot(), device="tiff", height=6.74, width=8.9, units="in", dpi=320, bg="white")

######## Background vs miRNA expression
library(DESeq2)
panelB=readRDS('/scratch/mmccall2_lab/ERCC_UMI/panel_B_SE.rds')
ddsColl=collapseReplicates(panelB, panelB$GEO_Sample)
raw_counts=assay(ddsColl)
total_counts=colSums(raw_counts)

cpm_sub=t(raw_counts)/(total_counts/1000000)
median_cpm_sub=apply(cpm_sub, 2, median)

source=rep("miRNA", length(median_cpm_sub))
idx=which(rowData(panelB)$miRNA=="-")
source[idx]="background (non-miRNA)"


df=data.frame(log(median_cpm_sub), source)
names(df)=c("log(median CPM)", "source")
ggplot(df, aes(x=`log(median CPM)`, fill=source))+geom_density(alpha=0.25)+ggtitle("Expression in sequences mapping to miRNA vs background expression")+
theme(plot.title=element_text(hjust=0.5, size=15), legend.position="bottom", legend.text=element_text(size=15),
      axis.title=element_text(size=15), axis.text=element_text(size=15), legend.title=element_text(size=15))

#ggsave("figures/supfigure3_B.tif", plot=last_plot(), device="tiff", height=6.74, width=8.9, units="in", dpi=320, bg="white")


ggarrange(p1, p2, ncol=2, nrow=1)
ggsave("figures/supfigure3.tif", plot=last_plot(), device="tiff", width=5.9, height=7.27, units="in", dpi=320, bg="white")

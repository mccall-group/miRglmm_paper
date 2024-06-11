## figure 2 for miRNA grant

## microRNAome data
library(microRNAome)
data(microRNAome)

## select just cell type data
icell <- which(microRNAome$sample_category == "cell_type")
dat <- microRNAome[ ,icell]

## filter samples with less than 100,000 total counts
i_rm <- which(colSums(assay(dat)) < 1e5)
dat <- dat[ ,-i_rm]

## proportion of non-zero counts
p_non_zero <- colSums(assay(dat)>0)

## proportion of reads going to each miRNA
p_reads <- sweep(assay(dat), 2, colSums(assay(dat)), FUN = "/")

## proportion of counts assigned to the top n miRNAs
prop_top <- function(d, output="prop"){
  d <- sort(d, decreasing = TRUE)
  if(output == "prop") return(cumsum(d / sum(d)))
  if(output == "name") return(names(d))
}
p_top_n <- apply(assay(dat), 2, prop_top)
top_n <- apply(assay(dat), 2, prop_top, "name")

## number of miRNAs taking up 90% of reads
n90 <- apply(p_top_n, 2, function(x) which(x > 0.9)[1])

## select cell types with at least 40 samples
library(stringr)
unique_celltypes <- names(which(table(dat$cell_tissue) > 40)) # 8 of them
ind <- which(dat$cell_tissue %in% unique_celltypes)
celltypes <- gsub("_", " ", dat$cell_tissue[ind], fixed = TRUE) 
celltypes <- str_replace(celltypes, "Endothelial", "Endothelial Cell")
celltypes <- str_replace(celltypes, "Smooth muscle", "Smooth Muscle Cell")
celltypes <- str_replace(celltypes, "Natural Killer Cells CD56", 
                         "Natural Killer Cell CD56")
celltypes <- str_replace(celltypes, "lymphocyte", 
                         "Lymphocyte")
celltypes <- factor(celltypes, levels = unique(celltypes))
df1 <- data.frame(celltype = celltypes,
                  n_mirna90 = n90[ind])

## stratify by cell type
library(ggplot2)
library(ggpubr)

fig2a <- ggplot(df1, aes(x=celltype, y=n_mirna90)) + geom_boxplot() +
  ylab("Number of miRNAs to which\n90% of reads are mapped") +
  xlab("Cell Type") + labs(tag = "A") +
  theme(plot.tag = element_text(size = 24, 
                                margin = margin(r = -20), 
                                face = "bold")) +
  theme(axis.title.x = element_text(size = 16)) + 
  theme(axis.title.y = element_text(size = 16)) + 
  theme(axis.text = element_text(size=14),
        axis.text.x = element_text(angle = 45, hjust = 1))

#############################

## mRNA data
library(SummarizedExperiment)
load("~/data/recount2_data/rseData_01_06_2020.RData")

## remove technical replicates from DRP001797
irm <- which(colData(rseData_01_06_2020)$project == "DRP001797" &
               duplicated(colData(rseData_01_06_2020)$sample))
rseData_01_06_2020 <- rseData_01_06_2020[, -irm]

## remove samples with average read length > 200
irm <- which(rseData_01_06_2020$avg_read_length > 200)
rseData_01_06_2020 <- rseData_01_06_2020[, -irm]

## total reads for each sample
nreads <- colSums(assay(rseData_01_06_2020))

## proportion of non-zero counts
p_non_zero <- colMeans(assay(rseData_01_06_2020) > 0)

## proportion of reads going to each mRNA
p_reads <- sweep(assay(rseData_01_06_2020), 2, 
                 colSums(assay(rseData_01_06_2020)), FUN = "/")

## proportion of counts assigned to the top n genes
prop_top <- function(d, output="prop"){
  d <- sort(d, decreasing = TRUE)
  if(output == "prop") return(cumsum(d))
  if(output == "name") return(names(d))
}
p_top_n <- matrix(nrow=nrow(p_reads), ncol=ncol(p_reads))
for(k in 1:ncol(p_reads)){
  p_top_n[ ,k] <- prop_top(p_reads[ ,k])
  if(k %% 100 == 0) cat(k/100)
}

## number of mRNAs taking up 90% of reads
n90 <- apply(p_top_n, 2, function(x) which(x > 0.9)[1])

## select celltypes with at least 10 samples
unique_celltypes <- names(which(table(rseData_01_06_2020$celltype) > 10)) # 7 of them
ind <- which(rseData_01_06_2020$celltype %in% unique_celltypes)
celltypes <- gsub("([a-z])([A-Z])", "\\1 \\2", 
                  rseData_01_06_2020$celltype[ind], 
                  perl = TRUE)
celltypes <- str_replace(celltypes, "Lymphocytes", 
                         "Lymphocyte")
df2 <- data.frame(celltype = celltypes,
                  n_mrna90 = n90[ind])

## stratify by cell type
library(ggplot2)
library(ggpubr)

fig2b <- ggplot(df2, aes(x=celltype, y=n_mrna90)) + geom_boxplot() +
  ylab("Number of mRNAs to which\n90% of reads are mapped") +
  xlab("Cell Type") + labs(tag = "B") +
  theme(plot.tag = element_text(size = 24, 
                                margin = margin(r = -20), 
                                face = "bold")) +
  theme(axis.title.x = element_text(size = 16)) + 
  theme(axis.title.y = element_text(size = 16)) + 
  theme(axis.text = element_text(size=14),
        axis.text.x = element_text(angle = 45, hjust = 1))

##################

fig2 <- ggarrange(fig2a, fig2b, heights = 1, widths = 1, 
                  nrow = 1, ncol = 2, align = "h")
#ggsave("~/Dropbox/grants/R01_microRNA/figure2.png", 
#       plot = fig2, device = "png", width = 12, height = 8)


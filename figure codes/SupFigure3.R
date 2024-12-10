library(pheatmap)
library(gridExtra)
library(reshape2)
# using supplementary data from https://pubmed.ncbi.nlm.nih.gov/33627690/ 

# Load miRNA-seq data, processed with sRNAbench from SRA study SRP090487
data <- read.table("miRNA-seq_RPM.tsv", header = TRUE, sep="\t", stringsAsFactors = FALSE)

# Set the first column as row names and remove it from the data
rownames(data) <- data[, 1]
data <- data[, -1]

# Calculate average expression and sort the rows
data$average_expression <- rowMeans(data, na.rm = TRUE)
sorted_data <- data[order(-data$average_expression), ]  # Sort in descending order

# Select the top 10 miRNAs based on average expression
top10_data <- sorted_data[1:10, -ncol(sorted_data)]  # Exclude the average_expression column
top5_data <- sorted_data[1:5, -ncol(sorted_data)]  # Exclude the average_expression column

# Define the vector of miRNAs of interest (top 5-10 most expressed)
#top 10
miRNAs_of_interest <- c("hsa-miR-486-5p", "hsa-miR-22-3p", "hsa-miR-92a-3p", "hsa-miR-16-5p", 
                        "hsa-miR-26a-5p", "hsa-miR-191-5p", "hsa-miR-126-5p", 
                        "hsa-miR-423-5p", "hsa-miR-451a", "hsa-miR-30d-5p")
#top 5
miRNAs_of_interest <- c("hsa-miR-486-5p", "hsa-miR-22-3p", "hsa-miR-92a-3p", "hsa-miR-16-5p", 
                        "hsa-miR-26a-5p")

# Calculate correlation among the top 10 miRNAs, miRNA to miRNA across samples
correlation_matrix_miRNAs <- cor(t(top10_data))



# Create the miRNA-seq heatmap and store it as a grob
heatmap1 <- pheatmap(correlation_matrix_miRNAs,
                     color = colorRampPalette(c("blue", "white", "red"))(100),
                     breaks = seq(-1, 1, length.out = 101),
                     border_color = NA,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     silent = TRUE,
                     main = "miRNA-seq")


# plot all qPCR + miRNA-seq

create_correlation_heatmap <- function(data, miRNAs_of_interest, title="18") {
  # Subset the data to include only the specified miRNAs
  subset_data <- data[rownames(data) %in% miRNAs_of_interest, ]
  # Reorder the correlation matrix according to the specified order
  ordered_data <- subset_data[miRNAs_of_interest, ]
  
  # Calculate correlation
  correlation_matrix <- cor(t(ordered_data))
  # print(correlation_matrix)
  
  # Generate heatmap and return as grob for further use
  heatmap <- pheatmap(correlation_matrix,
                      color = colorRampPalette(c("blue", "white", "red"))(100),
                      breaks = seq(-1, 1, length.out = 101),
                      border_color = NA,
                      cluster_rows = FALSE,
                      cluster_cols = FALSE,
                      silent = TRUE,
                      main = title)
  return(heatmap$gtable)
}

#top 10
miRNAs_of_interest <- c("hsa-miR-486-5p", "hsa-miR-22-3p", "hsa-miR-92a-3p", "hsa-miR-16-5p", 
                       "hsa-miR-26a-5p", "hsa-miR-191-5p", "hsa-miR-126-5p", 
                       "hsa-miR-423-5p", "hsa-miR-451a", "hsa-miR-30d-5p")
#top 5
miRNAs_of_interest <- c("hsa-miR-486-5p", "hsa-miR-22-3p", "hsa-miR-92a-3p", "hsa-miR-16-5p", 
                        "hsa-miR-26a-5p")

# List of qPCR datasets
qPCR_ABI <- read.delim2("ABI.txt", row.names=1)
qPCR_Exiqon <- read.delim2("Exiqon.txt", row.names=1)


datasets <- list(qPCR_ABI=qPCR_ABI, qPCR_Exiqon=qPCR_Exiqon)

heatmaps <- lapply(names(datasets), function(name) {
  create_correlation_heatmap(datasets[[name]], miRNAs_of_interest, name)
})

heatmaps[["miR"]] <- heatmap1$gtable

pdf("qPCR_Heatmaps.pdf", width = 16, height = 8)
grid.arrange(grobs = heatmaps, ncol = 2)
dev.off()
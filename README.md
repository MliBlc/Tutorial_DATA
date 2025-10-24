## Basic R scripts for begginers ##


# Load library
library(DESeq2)

# Load Count and meta data

count <- read.csv("RNA_seq_counts.txt" , header = T, sep = "\t", row.names = 1)

meta <- read.csv("metaTut.txt" , header = T, sep = "\t")

dds <- DESeqDataSetFromMatrix(count, meta, design = ~ Treatment)

smallestGroupSize <- 3

keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize

dds <- dds[keep,]

DESeq(dds) -> dds

# Remove the intermediate data

rm(count)

rm(meta)

rm(keep)

# Get All Results

Control_vs_Airpol_all <- results(object = dds, contrast = c("Treatment", "Airpol","Control"))

Airpol_vs_Antioxi_all <- results(object = dds, contrast = c("Treatment","Antioxi","Airpol"))


# Get up and down regulated genes

Control_vs_Airpol_up <- subset(Control_vs_Airpol_all, padj < 0.05 & log2FoldChange >= 1)

Control_vs_Airpol_down <- subset(Control_vs_Airpol_all, padj < 0.05 & log2FoldChange <= -1)

Airpol_vs_Antioxi_up <- subset(Airpol_vs_Antioxi_all, padj < 0.05 & log2FoldChange >= 1)

Airpol_vs_Antioxi_down <- subset(Airpol_vs_Antioxi_all, padj < 0.05 & log2FoldChange <= -1)

# Save results as csv
write.csv(Control_vs_Airpol_down, "Control_vs_Airpol_downregulated.csv")

write.csv(Control_vs_Airpol_up, "Control_vs_Airpol_upregulated.csv")

write.csv(Control_vs_Airpol_all, "Control_vs_Airpol_All_results.csv")

## Do the same for each comparison.

## Requirements
## 1. Expression count matrix derived from featurecounting algorithm after the alignment
## 2. Meta file contains grouping info of the samples as tab-seperated or comma seperated values.
## 3. R v4.1 or higher
## 4. R studio Latest version
## 5. DESeq2, EdgeR and LimmaR libraries




##Set the working directory


##Install BioManager and libraries
install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
BiocManager::install("limma")




#  $$$$$$$\  $$$$$$$$\  $$$$$$\  $$$$$$$$\  $$$$$$\   $$$$$$\  
#  $$  __$$\ $$  _____|$$  __$$\ $$  _____|$$  __$$\ $$  __$$\ 
#  $$ |  $$ |$$ |      $$ /  \__|$$ |      $$ /  $$ |\__/  $$ |
#  $$ |  $$ |$$$$$\    \$$$$$$\  $$$$$\    $$ |  $$ | $$$$$$  |
#  $$ |  $$ |$$  __|    \____$$\ $$  __|   $$ |  $$ |$$  ____/ 
#  $$ |  $$ |$$ |      $$\   $$ |$$ |      $$ $$\$$ |$$ |      
#  $$$$$$$  |$$$$$$$$\ \$$$$$$  |$$$$$$$$\ \$$$$$$ / $$$$$$$$\ 
#  \_______/ \________| \______/ \________| \___$$$\ \________|
#                                               \___|          
#                                                              
#                           


library(DESeq2)

# Load meta and count data
meta <- read.table("Meta.txt", header = TRUE, sep = "\t", stringsAsFactors = TRUE)
counts <- read.table("Counts.txt", header = TRUE, sep = "\t", row.names = 1)

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ treatments)


# Filter out genes with low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


# Run DESeq normalization and analysis
dds <- DESeq(dds)

# Set thresholds
fdr_threshold <- 0.05
logfc_threshold <- 1

### 1. Control vs Treatment
res1 <- results(dds, contrast = c("treatments", "Control", "Treatment"))
res1_df <- subset(res1)

# Subset for upregulated and downregulated genes
upregulated1 <- subset(res1, padj <= fdr_threshold & log2FoldChange >= logfc_threshold)
downregulated1 <- subset(res1, padj <= fdr_threshold & log2FoldChange <= -logfc_threshold)

# Save Control vs Treatment results
write.csv(res1_df, "Control_vs_Treatment_all_results.txt", quote = FALSE, row.names = TRUE)
write.csv(upregulated1, "Control_vs_Treatment_upregulated.txt", quote = FALSE, row.names = TRUE)
write.csv(downregulated1, "Control_vs_Treatment_downregulated.txt", quote = FALSE, row.names = TRUE)

### 2. Control vs Treatment2
res2 <- results(dds, contrast = c("treatments", "Control", "Treatment2"))
res2_df <- subset(res2)

# Subset for upregulated and downregulated genes

upregulated2 <- subset(res2, padj <= fdr_threshold & log2FoldChange >= logfc_threshold)
downregulated2 <- subset(res2, padj <= fdr_threshold & log2FoldChange <= -logfc_threshold)

# Save Control vs Treatment2 results
write.csv(res2_df, "Control_vs_Treatment2_all_results.txt", quote = FALSE, row.names = TRUE)
write.csv(upregulated2, "Control_vs_Treatment2_upregulated.txt", quote = FALSE, row.names = TRUE)
write.csv(downregulated2, "Control_vs_Treatment2_downregulated.txt", quote = FALSE, row.names = TRUE)

### 3. Treatment vs Treatment2
res3 <- results(dds, contrast = c("treatments", "Treatment2", "Treatment"))
res3_df <- subset(res3)

# Subset for upregulated and downregulated genes
upregulated3 <- subset(res3, padj <= fdr_threshold & log2FoldChange >= logfc_threshold)
downregulated3 <- subset(res3, padj <= fdr_threshold & log2FoldChange <= -logfc_threshold)

# Save Treatment vs Treatment2 results
write.csv(res3_df, "Treatment_vs_Treatment2_all_results.txt", quote = FALSE, row.names = TRUE)
write.csv(upregulated3, "Treatment_vs_Treatment2_upregulated.txt", quote = FALSE, row.names = TRUE)
write.csv(downregulated3, "Treatment_vs_Treatment2_downregulated.txt", quote = FALSE, row.names = TRUE)


#   .----------------.  .----------------.  .----------------.  .----------------.  .----------------. 
#  | .--------------. || .--------------. || .--------------. || .--------------. || .--------------. |
#  | |  _________   | || |  ________    | || |    ______    | || |  _________   | || |  _______     | |
#  | | |_   ___  |  | || | |_   ___ `.  | || |  .' ___  |   | || | |_   ___  |  | || | |_   __ \    | |
#  | |   | |_  \_|  | || |   | |   `. \ | || | / .'   \_|   | || |   | |_  \_|  | || |   | |__) |   | |
#  | |   |  _|  _   | || |   | |    | | | || | | |    ____  | || |   |  _|  _   | || |   |  __ /    | |
#  | |  _| |___/ |  | || |  _| |___.' / | || | \ `.___]  _| | || |  _| |___/ |  | || |  _| |  \ \_  | |
#  | | |_________|  | || | |________.'  | || |  `._____.'   | || | |_________|  | || | |____| |___| | |
#  | |              | || |              | || |              | || |              | || |              | |
#  | '--------------' || '--------------' || '--------------' || '--------------' || '--------------' |
#   '----------------'  '----------------'  '----------------'  '----------------'  '----------------' 


library(edgeR)

# Load meta and count data
meta <- read.table("Meta.txt", header = TRUE, sep = "\t", row.names = 1)
counts <- read.table("Counts.txt", header = TRUE, sep = "\t", row.names = 1)

# Create DGEList object
group <- factor(meta$treatments)
dge <- DGEList(counts = counts, group = group)

# Filter lowly expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Normalize the data
dge <- calcNormFactors(dge)

# Design matrix for the comparisons
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Estimate dispersion
dge <- estimateDisp(dge, design)

# Fit the model
fit <- glmFit(dge, design)

# Set thresholds
fdr_threshold <- 0.05
logfc_threshold <- 1

### 1. Control vs Treatment
contrast1 <- makeContrasts(Treatment_vs_Control = Treatment - Control, levels = design)
lrt1 <- glmLRT(fit, contrast = contrast1)
results1 <- topTags(lrt1, n=Inf)$table

upregulated1 <- results1[results1$logFC > logfc_threshold & results1$FDR <= fdr_threshold, ]
downregulated1 <- results1[results1$logFC < -logfc_threshold & results1$FDR <= fdr_threshold, ]

# Save Control vs Treatment results
write.table(results1, "Control_vs_Treatment_all_results.txt", sep = "\t", quote = FALSE, row.names = T)
write.table(upregulated1, "Control_vs_Treatment_upregulated.txt", sep = "\t", quote = FALSE, row.names = T)
write.table(downregulated1, "Control_vs_Treatment_downregulated.txt", sep = "\t", quote = FALSE, row.names = T)

### 2. Control vs Treatment2
contrast2 <- makeContrasts(Treatment2_vs_Control = Treatment2 - Control, levels = design)
lrt2 <- glmLRT(fit, contrast = contrast2)
results2 <- topTags(lrt2, n=Inf)$table

upregulated2 <- results2[results2$logFC > logfc_threshold & results2$FDR <= fdr_threshold, ]
downregulated2 <- results2[results2$logFC < -logfc_threshold & results2$FDR <= fdr_threshold, ]

# Save Control vs Treatment2 results
write.table(results2, "Control_vs_Treatment2_all_results.txt", sep = "\t", quote = FALSE, row.names = T)
write.table(upregulated2, "Control_vs_Treatment2_upregulated.txt", sep = "\t", quote = FALSE, row.names = T)
write.table(downregulated2, "Control_vs_Treatment2_downregulated.txt", sep = "\t", quote = FALSE, row.names = T)

### 3. Treatment vs Treatment2
contrast3 <- makeContrasts(Treatment2_vs_Treatment = Treatment2 - Treatment, levels = design)
lrt3 <- glmLRT(fit, contrast = contrast3)
results3 <- topTags(lrt3, n=Inf)$table

upregulated3 <- results3[results3$logFC > logfc_threshold & results3$FDR <= fdr_threshold, ]
downregulated3 <- results3[results3$logFC < -logfc_threshold & results3$FDR <= fdr_threshold, ]

# Save Treatment vs Treatment2 results
write.table(results3, "Treatment_vs_Treatment2_all_results.txt", sep = "\t", quote = FALSE, row.names = T)
write.table(upregulated3, "Treatment_vs_Treatment2_upregulated.txt", sep = "\t", quote = FALSE, row.names = T)
write.table(downregulated3, "Treatment_vs_Treatment2_downregulated.txt", sep = "\t", quote = FALSE, row.names = T)


#   _      _                              
#  | |    (_)                             
#  | |     _  _ __ ___   _ __ ___    __ _ 
#  | |    | || '_ ` _ \ | '_ ` _ \  / _` |
#  | |____| || | | | | || | | | | || (_| |
#  \_____/|_||_| |_| |_||_| |_| |_| \__,_|
#                                         
#  


library(limma)

# Load meta and count data
meta <- read.table("Meta.txt", header = TRUE, sep = "\t", stringsAsFactors = TRUE)
counts <- read.table("Counts.txt", header = TRUE, sep = "\t", row.names = 1)

# Prepare design matrix for the comparisons
group <- factor(meta$treatments)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Apply voom transformation
v <- voom(counts, design)

# Fit the linear model
fit <- lmFit(v, design)

# Create contrasts for each comparison
contrast_matrix <- makeContrasts(
  Treatment_vs_Control = Treatment - Control,
  Treatment2_vs_Control = Treatment2 - Control,
  Treatment2_vs_Treatment = Treatment2 - Treatment,
  levels = design
)

# Apply contrasts to the fit object
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Set thresholds
fdr_threshold <- 0.05
logfc_threshold <- 1

### 1. Control vs Treatment
res1 <- topTable(fit2, coef = "Treatment_vs_Control", number = Inf)
upregulated1 <- res1[res1$logFC > logfc_threshold & res1$adj.P.Val <= fdr_threshold, ]
downregulated1 <- res1[res1$logFC < -logfc_threshold & res1$adj.P.Val <= fdr_threshold, ]

# Save Control vs Treatment results
write.table(res1, "Control_vs_Treatment_all_results.txt", sep = "\t", quote = FALSE, row.names = T)
write.table(upregulated1, "Control_vs_Treatment_upregulated.txt", sep = "\t", quote = FALSE, row.names = T)
write.table(downregulated1, "Control_vs_Treatment_downregulated.txt", sep = "\t", quote = FALSE, row.names = T)

### 2. Control vs Treatment2
res2 <- topTable(fit2, coef = "Treatment2_vs_Control", number = Inf)
upregulated2 <- res2[res2$logFC > logfc_threshold & res2$adj.P.Val <= fdr_threshold, ]
downregulated2 <- res2[res2$logFC < -logfc_threshold & res2$adj.P.Val <= fdr_threshold, ]

# Save Control vs Treatment2 results
write.table(res2, "Control_vs_Treatment2_all_results.txt", sep = "\t", quote = FALSE, row.names = T)
write.table(upregulated2, "Control_vs_Treatment2_upregulated.txt", sep = "\t", quote = FALSE, row.names = T)
write.table(downregulated2, "Control_vs_Treatment2_downregulated.txt", sep = "\t", quote = FALSE, row.names = T)

### 3. Treatment vs Treatment2
res3 <- topTable(fit2, coef = "Treatment2_vs_Treatment", number = Inf)
upregulated3 <- res3[res3$logFC > logfc_threshold & res3$adj.P.Val <= fdr_threshold, ]
downregulated3 <- res3[res3$logFC < -logfc_threshold & res3$adj.P.Val <= fdr_threshold, ]

# Save Treatment vs Treatment2 results
write.table(res3, "Treatment_vs_Treatment2_all_results.txt", sep = "\t", quote = FALSE, row.names = T)
write.table(upregulated3, "Treatment_vs_Treatment2_upregulated.txt", sep = "\t", quote = FALSE, row.names = T)
write.table(downregulated3, "Treatment_vs_Treatment2_downregulated.txt", sep = "\t", quote = FALSE, row.names = T)

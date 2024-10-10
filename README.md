######## Prep ########
1. Set Working Directory
2. 
Use setwd() command to set the working directory where your data files are located.

3. Install Required Libraries
   
Use install.packages("BiocManager") command to install the BiocManager package, which helps in managing Bioconductor packages.

Use BiocManager::install("DESeq2") command to install the DESeq2 package for differential expression analysis.

Use BiocManager::install("edgeR") command to install the edgeR package for count data analysis.

Use BiocManager::install("limma") command to install the limma package for linear modeling and analysis.



4. Load Metadata and Count Data


Use read.table() command to load your metadata file (e.g., Meta.txt) into R.

Use read.table() command to load your count data file (e.g., Counts.txt) into R, specifying the correct delimiter and setting the row names appropriately.



####### DESeq2 Analysis ######

1a. Create a DESeqDataSet:

Use DESeqDataSetFromMatrix() command to create a DESeqDataSet from the count data and metadata.


2a. Filter Out Low Counts:

Use rowSums() command to filter out low-count genes based on a defined threshold.


3a. Run DESeq Analysis:

Use DESeq() command to run the DESeq analysis on the DESeqDataSet.


4a. Obtain Results for Comparisons:

Use results() command to extract results for specific comparisons.


5a. Export DESeq2 Results


Use write.table() or write.csv() command to save the results of the DEseq analysis.


######## edgeR Analysis ########


1b. Create a Factor for Treatment Groups:

Use factor() command to create a treatment group factor based on your metadata.


2b. Create a DGEList Object:

Use DGEList() command to create a DGEList object from the count data and treatment factor.


3b. Filter Lowly Expressed Genes:

Use filterByExpr() command to filter out lowly expressed genes.


4b. Normalize the Data:

Use calcNormFactors() command to normalize the data in your DGEList object.


5b. Create a Design Matrix:

Use model.matrix() command to create a design matrix based on treatment groups.


6b. Estimate Dispersion:

Use estimateDisp() command to estimate dispersion for the DGEList object.


7b. Fit the Model:
Use glmFit() command to fit the model to your data.


8b. Perform Likelihood Ratio Tests:
Use glmLRT() command to perform likelihood ratio tests for each comparison.


9b. Export edgeR Results
Use write.table() or write.csv() command to save the results of the EdgeR analysis.


################## limma Analysis #################

1c. Create a Design Matrix:
Define treatments as factors with factor(). Use model.matrix() command to create a design model matrix based on your treatment groups. Use colnames() to add treatment names as column names to design matrix.


group <- factor(meta$treatments)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)


2c. Apply the Voom Transformation:
Use voom() command to apply the voom transformation to your count data.


3c. Fit the Linear Model:
Use lmFit() command to fit the linear model to the transformed data.


4c. Create Contrasts for Each Comparison:
Use makeContrasts() command to define contrasts for your comparisons.


5c. Apply Empirical Bayes Moderation:
Use eBayes() command to apply empirical Bayes moderation to the fitted model.

6c. Obtain Results for Comparisons:
Use topTable() command to summarize the results for each contrast.

7c. Export limma Results
Use write.table() or write.csv() command to save the results of the limma analysis.

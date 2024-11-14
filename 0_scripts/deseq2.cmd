
Install R and DESeq2. Upon installing R, install DESeq2 on R:
> source("https://bioconductor.org/biocLite.R")
> biocLite("DESeq2")
Import DESeq2 library in R
> library("DESeq2")
Load gene(/transcript) count matrix and labels
> countData <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
> colData <- read.csv(PHENO_DATA, sep="\t", row.names=1)
Note: The PHENO_DATA file contains information on each sample, e.g., sex or population. The exact way to import this depends on the format of the file.
Check all sample IDs in colData are also in CountData and match their orders
> all(rownames(colData) %in% colnames(countData))
[1] TRUE
> countData <- countData[, rownames(colData)]
> all(rownames(colData) == colnames(countData))
[1] TRUE
Create a DESeqDataSet from count matrix and labels
> dds <- DESeqDataSetFromMatrix(countData = countData,
        colData = colData, design = ~ CHOOSE_FEATURE)
Run the default analysis for DESeq2 and generate results table
> dds <- DESeq(dds)
> res <- results(dds)
Sort by adjusted p-value and display
> (resOrdered <- res[order(res$padj), ]) 

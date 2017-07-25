
#SERIOUSLY Thanks to https://harshinamdar.wordpress.com/2014/11/11/quick-tutorial-on-deseq2/

library("DESeq2")

args <- commandArgs(trailingOnly = TRUE)

expr_input   = args[1];
sample_info  = args[2];
output_table = args[3];


countData <- read.table(expr_input,header=TRUE,row.names=1, sep="\t")
countData <- countData[grep("^__", row.names(countData), invert=TRUE),]
countData <- round(countData, 0)
colData <- read.table(sample_info,header=TRUE,row.names=1, sep="\t")

dds <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~ condition)
dds <- estimateSizeFactors(dds)
deseq_Ncounts <- counts(dds, normalized=TRUE)

write.table(deseq_Ncounts, output_table, sep="\t")

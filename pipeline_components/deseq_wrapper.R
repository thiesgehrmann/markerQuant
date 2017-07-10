#SERIOUSLY Thanks to https://harshinamdar.wordpress.com/2014/11/11/quick-tutorial-on-deseq2/

library("DESeq2")

args <- commandArgs(trailingOnly = TRUE)

expr_input   = args[1];
sample_info  = args[2];
conditionA   = args[3];
conditionB   = args[4];
output_table = args[5];

countData <- read.table(expr_input,header=TRUE,row.names=1, sep="\t")
# Remove the additional stuff from the marker and alignment quantification stuff (multimapped, etc)
countData <- countData[grep("^__", row.names(countData), invert=TRUE),]
countData <- round(countData, 0)

colData <- read.table(sample_info,header=TRUE,row.names=1, sep="\t")

dds <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~ condition)
dds <- DESeq(dds)
result_A_B <- results(dds, contrast=c("condition",conditionA, conditionB)) ## contrast specifies conditions to be tested
result_A_B <- result_A_B[complete.cases(result_A_B),] ## to remove rows with NA

write.table(result_A_B, output_table, sep="\t")






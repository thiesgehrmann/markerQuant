args <- commandArgs(trailingOnly = TRUE)

tests_input  = args[1];
tests_output = args[2];

tests <- read.table(tests_input,header=TRUE, sep="\t")
tests$padj <- p.adjust(tests$pvalue, method="BH")

write.table(tests, tests_output, row.names=FALSE, sep="\t")



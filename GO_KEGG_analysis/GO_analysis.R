# BiocManager::install("clusterProfiler")
rm(list=ls())
library("clusterProfiler")

term2gene <- read.csv("Zmays.GO.anno.term_gene.txt", header=F, sep="\t")
term2name <- read.csv("Zmays.GO.anno.term_name.txt", header=F, sep="\t")


gene <- read.csv("test_list.txt", header = F, sep="\t")
gene <- as.factor(gene$V1)
x <- enricher(gene, TERM2GENE=term2gene, TERM2NAME=term2name,
              pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)
write.table(x, "test_list.annotation.GO.txt", quote=F, row.names=F, col.names=T, sep="\t")


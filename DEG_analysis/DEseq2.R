##############################################################
# R code for PCA plot by DEseq2
#
# Xianqing Jia @Northwest University, jiaxq.nju@gmail.com
#
# 2023-10-09
##############################################################
#########################################################################
rm(list=ls())
wkdir           <- commandArgs(TRUE)[1]
input           <- commandArgs(TRUE)[2]
symlist         <- commandArgs(TRUE)[3]
replicate       <- commandArgs(TRUE)[4]
output          <- commandArgs(TRUE)[5]
#########################################################################

library(ggplot2)
library(ggsci)
library(DESeq2)
library(ggrepel)

#####################################################################################
## DEG analysis
setwd(wkdir)
mycounts<-read.table(input, header = T, row.names=1, sep = "\t")
mydata_all <- round(as.matrix(mycounts))

sym <- as.vector(unlist(strsplit(symlist, split = ",")))

mydata <- mydata_all[, sym]

condition <- factor(rep(c("T1", "T2"), each = replicate), levels = c("T1", "T2"))
mycolData <- data.frame(row.names=colnames(mydata), condition)

dds <- DESeqDataSetFromMatrix(countData = mydata, colData = mycolData, design =~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
res <- results(dds)

# table(res$padj <0.05)
# summary(res)

res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
write.table(resdata, file = output, sep="\t", quote=F, row.names=F, col.names=T)

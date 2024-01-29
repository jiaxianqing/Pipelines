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
size            <- commandArgs(TRUE)[5]
output          <- commandArgs(TRUE)[6]
#########################################################################

library(ggplot2)
library(ggsci)
library(DESeq2)
library(ggrepel)

#####################################################################################
## PCA Plot
setwd(wkdir)
mycounts<-read.table(input, header = T, row.names=1, sep = "\t")
mydata <- round(as.matrix(mycounts))

sym <- as.vector(unlist(strsplit(symlist, split = ",")))
isize <- as.vector(unlist(strsplit(size,split = ",")))

condition <- factor(rep(sym, each = replicate), levels = sym)
mycolData <- data.frame(row.names=colnames(mydata), condition)

dds <- DESeqDataSetFromMatrix(countData = mydata, colData = mycolData, design =~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)

vsd <- vst(dds)

pdf(file=output, width = as.numeric(isize[1]), height = as.numeric(isize[2]))

plotPCA(vsd) + 
  theme_bw() + 
  #scale_color_futurama() +
  geom_text_repel(aes(label = group), size = 3)

dev.off()

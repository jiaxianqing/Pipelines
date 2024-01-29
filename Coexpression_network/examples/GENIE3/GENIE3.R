# https://www.jianshu.com/p/d71dcd4cff5a
# https://bioconductor.riken.jp/packages/release/bioc/vignettes/GENIE3/inst/doc/GENIE3.html

# BiocManager::install("GENIE3")
# install.packages("foreach")
# install.packages("doParallel")
# install.packages("doRNG")
library(GENIE3)


set.seed(123)

rm(list = ls())
fpkm <- read.table('L_all.count.DEseq2.xls.DEG.padj001.list.count.xls', header = T, row.names = 1)
list <- read.table('L.TF.list.txt', header = T)

mydata <- round(as.matrix(fpkm))

weightList <- GENIE3(mydata, nCores=4, targets=names(list$gene_id), regulators=list$gene_id)


# TFlinkList <- getLinkList(weightList, threshold=0.015)

TFlinkList <- getLinkList(weightList, reportMax=1000)
write.table(TFlinkList,"L_all.count.DEseq2.xls.DEG.padj001.list.count.GENIE3.TF_top1000.xls",
            quote=F,row.names=F,col.names=T,sep="\t")
TFlinkList <- getLinkList(weightList, reportMax=2000)
write.table(TFlinkList,"L_all.count.DEseq2.xls.DEG.padj001.list.count.GENIE3.TF_top2000.xls",
            quote=F,row.names=F,col.names=T,sep="\t")

###
# Shoot
rm(list = ls())
fpkm <- read.table('S_all.count.DEseq2.xls.DEG.padj001.list.count.xls', header = T, row.names = 1)
list <- read.table('S.TF.list.txt', header = T)

mydata <- round(as.matrix(fpkm))

weightList <- GENIE3(mydata, nCores=4, targets=names(list$gene_id), regulators=list$gene_id)

TFlinkList <- getLinkList(weightList, reportMax=1000)
write.table(TFlinkList,"S_all.count.DEseq2.xls.DEG.padj001.list.count.GENIE3.TF_top1000.xls",
            quote=F,row.names=F,col.names=T,sep="\t")
TFlinkList <- getLinkList(weightList, reportMax=2000)
write.table(TFlinkList,"S_all.count.DEseq2.xls.DEG.padj001.list.count.GENIE3.TF_top2000.xls",
            quote=F,row.names=F,col.names=T,sep="\t")


###
# Root
rm(list = ls())
fpkm <- read.table('R_all.count.DEseq2.xls.DEG.padj001.list.count.xls', header = T, row.names = 1)
list <- read.table('R.TF.list.txt', header = T)

mydata <- round(as.matrix(fpkm))

weightList <- GENIE3(mydata, nCores=4, targets=names(list$gene_id), regulators=list$gene_id)

TFlinkList <- getLinkList(weightList, reportMax=1000)
write.table(TFlinkList, "R_all.count.DEseq2.xls.DEG.padj001.list.count.GENIE3.TF_top1000.xls",
            quote=F, row.names=F, col.names=T, sep="\t")
TFlinkList <- getLinkList(weightList, reportMax=2000)
write.table(TFlinkList, "R_all.count.DEseq2.xls.DEG.padj001.list.count.GENIE3.TF_top2000.xls",
            quote=F, row.names=F, col.names=T, sep="\t")
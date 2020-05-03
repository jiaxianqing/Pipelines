##############################################################
# R code for DEG analysis
#
# https://www.nature.com/articles/nprot.2016.095
#
# 2018-1-25
##############################################################
#########################################################################
wkdir           <- commandArgs(TRUE)[1]
sample_list     <- commandArgs(TRUE)[2]
pheno_file      <- commandArgs(TRUE)[3]
outtrans        <- commandArgs(TRUE)[4]
outgene         <- commandArgs(TRUE)[5]
#########################################################################
#source("http://bioconductor.org/biocLite.R")
#biocLite("ballgown")
#biocLite("RSkittleBrewer")
#biocLite("devtools")
#devtools::install_github('alyssafrazee/RSkittleBrewer', force = TRUE)

library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

setwd(wkdir)
list <- read.table(sample_list, header = F)
sampl <- as.vector(list$V1)
bgdata <- ballgown(samples=sampl)

pheno_data <- read.delim(pheno_file)
pData(bgdata) = pheno_data[match(sampleNames(bgdata), pheno_data$sampleID), ]

bgdata.filt <- subset(bgdata, "rowVars(texpr(bgdata)) > 1", genomesubset=TRUE)

results_transcripts = stattest(bgdata.filt, feature="transcript",covariate="group", getFC=TRUE, meas="FPKM")
results_transcripts = data.frame(geneNames=ballgown::geneNames(bgdata.filt), geneIDs=ballgown::geneIDs(bgdata.filt), results_transcripts)
results_genes = stattest(bgdata.filt, feature="gene", covariate="group", getFC=TRUE, meas="FPKM")
results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes, pval)

write.table(results_transcripts, outtrans, sep = "\t",quote = F, row.names=FALSE)
write.table(results_genes, outgene, sep = "\t",quote = F, row.names=FALSE)


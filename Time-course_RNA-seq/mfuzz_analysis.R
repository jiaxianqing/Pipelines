##
# BiocManager::install("Mfuzz")
library("Mfuzz")
rm(list=ls())

tdata <- read.table("all_deg.fc1_5.T1_T2.fpkm.avr.xls", header = T, row.names= 1, sep = "\t")

gene_tpm <- data.matrix(tdata)
eset <- new("ExpressionSet",exprs = gene_tpm)

gene.r <- filter.NA(eset, thres=0.25)
gene.f <- fill.NA(gene.r,mode="mean")
# gene.f <- fill.NA(gene.r,mode="knn")

tmp <- filter.std(gene.f,min.std=0)
gene.s <- standardise(tmp)

test = as.data.frame(gene.s@assayData[["exprs"]])
write.table(test,"all_deg.fc1_5.T1_T2.fpkm.avr.expression_changes.txt",quote=F,row.names=T,col.names=T,sep="\t")



### test the optimal cluster number 
## method1 -> 15
m <- mestimate(gene.s)
# tmp <- cselection(gene.s, m=m, crange=seq(2,40,2), repeats=3, visu=TRUE)
Dmin(gene.s, m, crange=seq(2,20,1), repeats=3, visu=TRUE)

### method2 -> 12
# install.packages('vegan', dependencies = TRUE)
library(vegan)
ca_clust <- cascadeKM(tdata, 1, 20, iter = 1000)
ca_clust$results

calinski.best <- as.numeric(which.max(ca_clust$results[2,]))
calinski.best
pdf("vegan.k20.pdf")
plot(ca_clust, sortg = TRUE, grpmts.plot = TRUE)
dev.off()


# 
c <- 12
m <- mestimate(gene.s)

cl <- mfuzz(gene.s, c = c, m = m)
cl$size

# cl$cluster[cl$cluster == 1]
pdf("all_deg.fc1_5.T1_T2.fpkm.avr.k17.pdf",width = 10, height=6)
mfuzz.plot(gene.s,cl,mfrow=c(2,3), new.window= FALSE)
dev.off()

    # mfuzz.plot2(gene.s,cl,mfrow=c(1,1),min.mem=0,
#             ylim.set=c(0,0), xlab="Time",ylab="Expression changes",x11=TRUE,
#             ax.col="black",bg = "white",col.axis="black",col.lab="black",
#             col.main="black",col.sub="black",col="black",centre=FALSE,
#             centre.col="black",centre.lwd=10)

sig_gene_info = as.data.frame(cl$cluster)
write.table(sig_gene_info,"all_deg.fc1_5.T1_T2.fpkm.avr.k12.cluster_info.txt",quote=F,row.names=T,col.names=F,sep="\t")



O <- overlap(cl)
Ptmp <-  overlap.plot(cl,over=O,thres=0.05)
cor(t(cl[[1]]))
acore <- acore(gene.s,cl,min.acore=0)
acore_list <- do.call(rbind, lapply(seq_along(acore), function(i){ data.frame(CLUSTER=i, acore[[i]])}))

write.table(acore_list,"acore_list.output.C6.txt",quote=F,row.names=F,col.names=T,sep="\t")

## Time-course RNA-seq   

> Dr. Xianqing Jia (jiaxianqing@nwu.edu.cn)   
> Latest update: 2024-01-29   
> Lab website: https://jialab.life/   

---
This pipeline is based on **Mfuzz** in R:   
http://mfuzz.sysbiolab.eu/   
https://www.bioconductor.org/packages/release/bioc/html/Mfuzz.html

```
BiocManager::install("Mfuzz")
```

It includes:
+ [Step 1. Data pre-processing](#step-1-data-pre-processing)   
+ [Step 2. Testing the optimal cluster number](#step-2-testing-the-optimal-cluster-number)   
+ [Step 3. Soft clustering](#step-3-soft-clustering)   
+ [Step 4. Global clustering structures](#step-4-global-clustering-structures)   
+ [Step 5. Visualization](#step-5-visualization)   

---
### Step 1. Data pre-processing
Usually, normalized expression profiles are used for Mfuzz anlysis, such as TMP values obtained by [StringTie](https://github.com/jiaxianqing/Pipelines/tree/master/DEG_analysis/#1-stringtie-and-ballgown).

```
library("Mfuzz")
rm(list=ls())

tdata <- read.table("all_deg.fc1_5.T1_T2.tpm.avr.xls", header = T, row.names= 1, sep = "\t")

gene_tpm <- data.matrix(tdata)
eset <- new("ExpressionSet",exprs = gene_tpm)

gene.r <- filter.NA(eset, thres=0.25)
gene.f <- fill.NA(gene.r, mode="mean")
# gene.f <- fill.NA(gene.r, mode="knn")

tmp <- filter.std(gene.f, min.std=0)
gene.s <- standardise(tmp)

test = as.data.frame(gene.s@assayData[["exprs"]])
write.table(test, "all_deg.fc1_5.T1_T2.tpm.avr.expression_changes.txt", quote=F,row.names=T, col.names=T, sep="\t")
```

### Step 2. Testing the optimal cluster number 
- Method1
```
m <- mestimate(gene.s)
# tmp <- cselection(gene.s, m=m, crange=seq(2,40,2), repeats=3, visu=TRUE)
Dmin(gene.s, m, crange=seq(2,20,1), repeats=3, visu=TRUE)
```

- Method2 (recommended)
```
#install.packages('vegan', dependencies = TRUE)
library(vegan)
ca_clust <- cascadeKM(tdata, 1, 20, iter = 1000)
ca_clust$results

calinski.best <- as.numeric(which.max(ca_clust$results[2,]))
calinski.best
pdf("vegan.k20.pdf")
plot(ca_clust, sortg = TRUE, grpmts.plot = TRUE)
dev.off()
```

### Step 3. Soft clustering 

```
c <- 12 # optimal cluster number
m <- mestimate(gene.s)

cl <- mfuzz(gene.s, c = c, m = m)
cl$size

# cl$cluster[cl$cluster == 1]
pdf("all_deg.fc1_5.T1_T2.tpm.avr.k17.pdf", width = 10, height=6)
mfuzz.plot(gene.s, cl, mfrow=c(2,3), new.window= FALSE)
dev.off()

sig_gene_info = as.data.frame(cl$cluster)
write.table(sig_gene_info,"all_deg.fc1_5.T1_T2.tpm.avr.k12.cluster_info.txt",quote=F,row.names=T,col.names=F,sep="\t")
```
### Step 4. Global clustering structures
An interesting feature of soft clustering is the overlap or coupling between clusters.
```
O <- overlap(cl)
Ptmp <-  overlap.plot(cl,over=O, thres=0.05)

cl3 <- mfuzz(gene.s, c=10, m=1.25)
mfuzz.plot(gene.s, cl=cl3, mfrow=c(2,3))
O3 <- overlap(cl3)
overlap.plot(cl3,over=O3, P=Ptmp, thres=0.05)

cor(t(cl[[1]]))
acore <- acore(gene.s, cl, min.acore=0)
acore_list <- do.call(rbind, lapply(seq_along(acore), function(i){ data.frame(CLUSTER=i, acore[[i]])}))
write.table(acore_list, "acore_list.output.C6.txt", quote=F, row.names=F, col.names=T, sep="\t")
```

### Step 5. Visualization
```
library(ggplot2)
library(reshape2)
library(ggthemes)
library(ggsci)

rm(list=ls())
adata <- read.table("all_deg.fc1_5.T1_T2.tpm.avr.zscore.txt", header = T, sep = "\t")
fdata <- melt(adata, id.vars = c('ID','Cluster'), variable.name = 'Treat', value.name = 'Value')

write.table(fdata,"all_deg.fc1_5.T1_T2.tpm.avr.zscore.melt.txt",
            quote=F,row.names=F,col.names=T,sep="\t")

ffdata <- read.table("all_deg.fc1_5.T1_T2.tpm.avr.zscore.melt.mod.txt", header = T, sep = "\t")

ggplot(ffdata) +
  geom_line(aes(x = Treat, y = Value, group =ID), color= "grey") +
  facet_grid(Cluster~Group, scales = "free", space="free")+
  geom_line(data=aa, aes(x=Treat, y=Value, group = Cluster, color = Group), size = 0.8)+
  theme_few()+
  scale_color_lancet()+
  theme(legend.position="none",strip.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 45, hjust = 0.4, color = "black", vjust = 0.5),
        axis.text.y = element_text(color = "black"))

ggsave(filename = "all_deg.fc1_5.T1_T2.tpm.avr.zscore.melt.mod.pdf", heigh = 15, width =3.5, device = "pdf")
```


An example for final visualization:   
<img src="https://github.com/jiaxianqing/Pipelines/blob/master/Phylogeny_construction/species_tree.png" width = "25%" height = "25%" div align = "center" />

---

### Related publications
More details could also be found in our publications:

1. Zhang, Y., Wang, L., Guo, Z., Xu, L., Zhao, H., Zhao, P., Ma, C., Yi, K. and Jia, X. (2022) [Revealing the underlying molecular basis of phosphorus recycling in the green manure crop Astragalus sinicus](https://www.sciencedirect.com/science/article/pii/S0959652622005625). *Journal of Cleaner Production*, 341, 130924.
2. Wang, X., Wang, B., Song, Z., Zhao, L., Ruan, W., Gao, Y., Jia, X. and Yi, K. (2022) [A spatial–temporal understanding of gene regulatory networks and NtARF-mediated regulation of potassium accumulation in tobacco](https://link.springer.com/article/10.1007/s00425-021-03790-2). *Planta*, 255, 1–15. 


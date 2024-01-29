## Co-expression network analysis

> Dr. Xianqing Jia (jiaxianqing@nwu.edu.cn)   
> Latest update: 2024-01-29   
> Lab website: https://jialab.life/   

Two popular co-expression network analysis pipelines are included:
+ [1. GENIE3](#1-genie3)   
+ [2. ARACNe-AP](#2-aracne-ap) (Recommended)   

---
## 1. GENIE3
This pipeline is based on **GENIE3** in R:   
https://doi.org/10.1371/journal.pone.0012776   
https://bioconductor.org/packages/release/bioc/html/GENIE3.html
https://bioconductor.org/packages/release/bioc/vignettes/GENIE3/inst/doc/GENIE3.html

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GENIE3")
```

### Pipeline
```
library(GENIE3)
set.seed(123)
rm(list = ls())

fpkm <- read.table('L_all.count.DEseq2.xls.DEG.padj001.list.count.xls', header = T, row.names = 1)
list <- read.table('L.TF.list.txt', header = T)

mydata <- round(as.matrix(fpkm))

weightList <- GENIE3(mydata, nCores=4, targets=names(list$gene_id), regulators=list$gene_id)

TFlinkList <- getLinkList(weightList, reportMax=1000)
# TFlinkList <- getLinkList(weightList, threshold=0.015)
write.table(TFlinkList,"L_all.count.DEseq2.xls.DEG.padj001.list.count.GENIE3.TF_top1000.xls",
            quote=F,row.names=F,col.names=T,sep="\t")
TFlinkList <- getLinkList(weightList, reportMax=2000)
write.table(TFlinkList,"L_all.count.DEseq2.xls.DEG.padj001.list.count.GENIE3.TF_top2000.xls",
            quote=F,row.names=F,col.names=T,sep="\t")
```

## 2. ARACNe-AP
This pipeline is based on **ARACNe-AP**:   
https://doi.org/10.1093/bioinformatics/btw216   
https://github.com/califano-lab/ARACNe-AP

### Pipeline

```
mkdir -p /05_coexpression2/seeds
# Step 1
java -Xmx5G -jar /ARACNe-AP/dist/aracne.jar \
    -e /05_coexpression2/Root659_tpm.txt \
    -o /05_coexpression2/seeds \
    -t /05_coexpression2/Root659_TF.txt \
    -p 1E-8 -s 1 --calculateThreshold 

# Step 2
for i in {1..100}
do
    java -Xmx5G -jar /ARACNe-AP/dist/aracne.jar \
        -e /05_coexpression2/Root659_tpm.txt \
        -o /05_coexpression2/seeds \
        -t /05_coexpression2/Root659_TF.txt \
        -p 1E-8 -s $i
done

# Step 3
java -Xmx5G -jar /ARACNe-AP/dist/aracne.jar \
    -o /05_coexpression2/seeds \
    --consolidate

```

## Related publications
More details could also be found in our publications:

1. Wang, X., Wang, B., Song, Z., Zhao, L., Ruan, W., Gao, Y., Jia, X. and Yi, K. (2022) [A spatial–temporal understanding of gene regulatory networks and NtARF-mediated regulation of potassium accumulation in tobacco](https://link.springer.com/article/10.1007/s00425-021-03790-2). *Planta*, 255, 1–15. 

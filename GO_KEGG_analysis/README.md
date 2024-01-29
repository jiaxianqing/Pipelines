## GO & KEGG analysis

> Dr. Xianqing Jia (jiaxianqing@nwu.edu.cn)   
> Latest update: 2024-01-29   
> Lab website: https://jialab.life/   

**Before this instruction begins, I must to say that I hate GO and KEGG analysis very much. In my deep mind, I think GO/KEGG analysis is tasteless and it is abusive in making someone's paper look informative.**

This pipeline includes:
+ [GO analysis](#go-analysis)   
    + [1. AgriGO v2](#1-agrigo-v2)   
    + [2. clusterProfiler](#2-clusterprofiler)   
+ [KEGG analysis](#kegg-analysis)   

---
## GO analysis
### 1. AgriGO v2
[AgriGO](http://systemsbiology.cau.edu.cn/agriGOv2/) is a decent tool for gene ontology (GO) analysis. It is online and fast. Here is a Rscripts for ploting AgriGO results:
```
   library(ggplot2)
   library(ggthemes)

   go_rich <- read.table("DEGs.txt", sep = "\t", header = T)

   go_rich$Type <- factor(go_rich$Type,levels=c('Up_regulation','Down_regulation'))
   
   ggplot(go_rich, aes(-1*log10(pvalue),Term)) +
     facet_grid(Type+term_type~.,scales = "free",space="free")+
     geom_point(aes(size=queryitem,color=-1*log10(FDR))) +
     scale_colour_gradient(low="blue",high="red") +
     theme_bw()+
     labs(color=expression(-log[10]("FDR")),size="Gene count",
          x=expression(paste(-log[10], "(", italic(P), " ", value, ")", sep="")), 
          y="GO items")
   
   ggsave(filename = "DEGs.pdf", heigh = 5, width =6.5, device = "pdf")

   ## An example for GO results file (tab-separated):
   # GO_acc	term_type	Term	queryitem	querytotal	bgitem	bgtotal	pvalue	FDR	Type
   # GO:0006629	P	lipid metabolic process	11	108	528	24075	1.20E-07	1.80E-05	Up_regulation
   # GO:0008610	P	lipid biosynthetic process	5	108	180	24075	4.80E-05	0.0035	Up_regulation
```


### 2. clusterProfiler
This pipeline for DEG analysis is based on **clusterProfiler**.   
https://doi.org/10.1016/j.xinn.2021.100141   
https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html

```
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("clusterProfiler")
```

#### Pipeline
```
    rm(list=ls())
    library("clusterProfiler")

    term2gene <- read.csv("Zmays.GO.anno.term_gene.txt", header=F, sep="\t")
    term2name <- read.csv("Zmays.GO.anno.term_name.txt", header=F, sep="\t")

    gene <- read.csv("test_list.txt", header = F, sep="\t")
    gene <- as.factor(gene$V1)
    x <- enricher(gene, TERM2GENE=term2gene, TERM2NAME=term2name,
                pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)
    write.table(x, "test_list.annotation.GO.txt", quote=F, row.names=F, col.names=T, sep="\t")
```

## KEGG analysis
[KOBAS](http://kobas.cbi.pku.edu.cn) is the best tool in my mind.   
https://academic.oup.com/nar/article/49/W1/W317/6292104

Of course, you can also try **clusterProfiler**.

---
### Related publications
More details could also be found in our publications:
1. Zhang, Y., Wang, L., Guo, Z., Xu, L., Zhao, H., Zhao, P., Ma, C., Yi, K. and Jia, X. (2022) [Revealing the underlying molecular basis of phosphorus recycling in the green manure crop Astragalus sinicus](https://www.sciencedirect.com/science/article/pii/S0959652622005625). *Journal of Cleaner Production*, 341, 130924.
2. Wang, X., Wang, B., Song, Z., Zhao, L., Ruan, W., Gao, Y., Jia, X. and Yi, K. (2022) [A spatial–temporal understanding of gene regulatory networks and NtARF-mediated regulation of potassium accumulation in tobacco](https://link.springer.com/article/10.1007/s00425-021-03790-2). *Planta*, 255, 1–15. 
3. Jia, X., Yu, L., Tang, M., Tian, D., Yang, S., Zhang, X., & Traw, M. B. (2020). [Pleiotropic changes revealed by in situ recovery of the semi-dwarf gene sd1 in rice](http://www.sciencedirect.com/science/article/pii/S0176161720300298). *Journal of Plant Physiology*, 248, 153141.


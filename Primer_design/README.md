>Pipeline for batch primer design

## Backgroud
---

awk 'BEGIN{OFS="\t"} {print $1,$2-1000,$2+1000,$3;}' /home/jxq/Data/rice/project/head_date2/07_genes/primer/candidate.site.list | \
    ```getfasta -name -bed - -fo /home/jxq/Data/rice/project/head_date2/07_genes/primer/candidate.site.list.ex1k.fa```

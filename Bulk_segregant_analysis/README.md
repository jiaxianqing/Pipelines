## NGS Bulk Segregant Analysis
Bulk Segregant Analysis (BSA) is an elegant method to identify DNA markers tightly linked to the causal gene for a given phenotype. Following a cross between parental lines showing contrasting phenotypes, the resulting F2 progeny are scored for segregation of the phenotype. Two bulked DNA samples are generated from the progeny showing contrasting phenotypes, and DNA markers exhibiting differences between the two bulks are screened. Two typical literatures:
>[1.Takagi, H. et al. QTL-seq: rapid mapping of quantitative trait loci in rice by whole genome resequencing of DNA from two bulked populations. The Plant Journal 74, 174–183 (2013).](https://onlinelibrary.wiley.com/doi/abs/10.1111/tpj.12105)</br>

>[2.Wang, C. et al. Dissecting a heterotic gene through GradedPool-Seq mapping informs a rice-improvement strategy. Nat Commun 10, 1–12 (2019).](https://www.nature.com/articles/s41467-019-11017-y)</br>

This pipeline is based on a R package [QTLseqr](https://acsess.onlinelibrary.wiley.com/doi/full/10.3835/plantgenome2018.01.0006).

---
* **Prepare input file for QTLseqr**

Phenotype data (groups of samples) and site (marker) information are needed. There are two strategies for site information, alelle frequency (AF) and read depth (DP, also called allele depth, AD). Here an example for read depth DP strategy:

```
  bgzip -dc samples.ug.m1.vcf.gz | \
      awk 'BEGIN{OFS="\t"} !/^#/ {print $1,$2-1,$2}' | tabix -h -p vcf -R - \
          samples.ug.vcf.gz | \
              perl vcf2qtlseqr.pl --merge-method DP \
                  --minor-allele-frequency 0.05 --minor-total-depth 100 \
                  -i - -g bulk.group.list \
                      > samples.ug.m1.dp.tab

  # Note: 1. "samples.ug.m1.vcf.gz" and "samples.ug.vcf.gz" are from QTL-mapping pipeline.
  #       2. An example for "bulk.group.list" file (tab-separated):
  #           #Sample	HD
  #           DE55	Early
  #           DE2	Early
  #           DE58	Late
  #           DE56	Late
  #       3. An example for output file:
  #           CHROM   POS     REF     ALT     AD_REF.Early    AD_ALT.Early    GQ.Early        AD_REF.Late     AD_ALT.Late     GQ.Late
  #           chr01   1203    T       C       185     147     77      155     257     84
  #           chr01   1249    A       C       191     155     80      171     256     86
  #           chr01   1266    G       A       186     162     82      181     263     85
```

* **Run QTLseqr analysis**

```
  Rscript QTLseqr.analysis.dp.R \
      samples.ug.m1.dp.tab \
      Late Early \ # group names
      100000 \ # window size
      samples.ug.m1.dp # prefix for output files
```

* **Plot QTLseqr results**</br>
Two important indexes, ΔSNP-index and Gprime, are used to measure linked DNA markers. </br>

&emsp; 1. ΔSNP-index
```
  awk 'BEGIN{OFS="\t"} !/^CHROM/ {print $2,$3,$17}' samples.ug.m1.dp.QTLs.sites_info.csv | \
      perl smooth_window_count.pl -s - \
          -b IRGSP-1.0_genome.w100k.bed --methods average \
          > samples.ug.m1.dp.QTLseqr.w100k.deltaSNP.stats
  
  awk 'BEGIN{OFS="\t"} !/^CHROM/ {print $2,$3,$22}' samples.ug.m1.dp.QTLs.sites_info.csv | \
      perl smooth_window_count.pl -s - \
          -b IRGSP-1.0_genome.w100k.bed --methods average | cut -f 7 | sed 's/Stats/CI95/' \
          > samples.ug.m1.dp.QTLseqr.w100k.CI95.stats
  
  awk 'BEGIN{OFS="\t"} !/^CHROM/ {print $2,$3,$23}' samples.ug.m1.dp.QTLs.sites_info.csv | \
      perl smooth_window_count.pl -s - \
          -b IRGSP-1.0_genome.w100k.bed --methods average | cut -f 7 | sed 's/Stats/CI99/' \
          > samples.ug.m1.dp.QTLseqr.w100k.CI99.stats
  
  paste -d "\t" samples.ug.m1.dp.QTLseqr.w100k.deltaSNP.stats \
      samples.ug.m1.dp.QTLseqr.w100k.CI95.stats \
      samples.ug.m1.dp.QTLseqr.w100k.CI99.stats \
      > samples.ug.m1.dp.QTLseqr.w100k.deltaSNP.stats.csv

  Rscript plot_smooth_windows.deltaSNP.R \
      samples.ug.m1.dp.QTLseqr.w100k.deltaSNP.stats.csv \
      IRGSP-1.0_genome.centromere.bed \
      deltaSNP-index \
      bulk.QTLseqr.w100k.deltaSNP.stats \
      samples.ug.m1.dp.QTLseqr.w100k.deltaSNP.stats.pdf
```

&emsp; 2. Gprime
```
  awk 'BEGIN{OFS="\t"} !/^CHROM/ {print $2,$3,$25}' samples.ug.m1.dp.QTLs.sites_info.csv | \
      perl /home/jxq/Data/scripts/my_scripts/QTL_analysis/smooth_window_count.pl -s - \
          -b /home/jxq/Data/rice/ref/IRGSP-1.0_genome.w100k.bed --methods average \
          > samples.ug.m1.dp.QTLseqr.w100k.Gprime.stats
  
  awk 'BEGIN{OFS="\t"} !/^CHROM/ {print $2,$3,$28}' samples.ug.m1.dp.QTLs.sites_info.csv | \
      perl /home/jxq/Data/scripts/my_scripts/QTL_analysis/smooth_window_count.pl -s - \
          -b /home/jxq/Data/rice/ref/IRGSP-1.0_genome.w100k.bed --methods average | \
          cut -f 7 | sed 's/Stats/qvalue/' \
          > samples.ug.m1.dp.QTLseqr.w100k.qvalue.stats
  
  paste -d "\t" samples.ug.m1.dp.QTLseqr.w100k.Gprime.stats \
      samples.ug.m1.dp.QTLseqr.w100k.qvalue.stats \
      > samples.ug.m1.dp.QTLseqr.w100k.Gprime.stats.csv

  Rscript /home/jxq/Data/scripts/my_scripts/QTL_analysis/plot_smooth_windows.Gprime.R \
      samples.ug.m1.dp.QTLseqr.w100k.Gprime.stats.csv \
      /home/jxq/Data/rice/ref/centromere.bed \
      Gprime 29s_bulk.QTLseqr.w100k.Gprime.stats 0.01 \
      samples.ug.m1.dp.QTLseqr.w100k.Gprime.stats.pdf
```

* **Find related genes**
```
  sed 's/\"//g' samples.ug.m1.dp.QTLs.deltaSNP.csv | \
      awk -F "," 'BEGIN{OFS="\t"} !/CHROM/ {print $1,$3,$4}' | \
      bedtools intersect -a - -b all.mod.gff3.genes.csv -wa -wb \
      > samples.ug.m1.dp.QTLs.deltaSNP.gene.msu.bed

  sed 's/\"//g' samples.ug.m1.dp.QTLs.Gprime.csv | \
      awk -F "," 'BEGIN{OFS="\t"} !/CHROM/ {print $1,$3,$4}' | \
      bedtools intersect -a - -b all.mod.gff3.genes.csv -wa -wb \
      > samples.ug.m1.dp.QTLs.Gprime.gene.msu.bed
```



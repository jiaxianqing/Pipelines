## Differentially expressed gene (DEG) analysis

This pipeline for DEG analysis is followed [Pertea’s protocol](https://www.nature.com/articles/nprot.2016.095): HISAT2, StringTie and Ballgown.

---

### Software dependencies:

* [HISAT2](http://daehwankimlab.github.io/hisat2/download/)
* [gffread](https://github.com/gpertea/gffread)
* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Picard](https://broadinstitute.github.io/picard/)
* [samtools](http://samtools.sourceforge.net/)
* [StringTie](https://ccb.jhu.edu/software/stringtie/)
* R packages:  
    * ballgown
    * RSkittleBrewer
    * genefilter
    * dplyr
    * devtools
---

### Pipeline
* **Genome alignment (mapping) by HISAT2**

&emsp; 1. Prepare reference files for HISAT2
```
  mkdir -p /msu7
  gffread /MSU_osa1r7/all.mod.gff3 -T \
    -o /MSU_osa1r7/all.mod.gtf
  extract_exons.py MSU_osa1r7/all.mod.gtf \
      > /msu7/all.mod.exon
  extract_splice_sites.py MSU_osa1r7/all.mod.gtf \
      > msu7/all.mod.ss
  hisat2-build -p 8 IRGSP-1.0_genome.fasta \
      --ss /msu7/all.mod.ss \
      --exon /msu7/all.mod.exon\
      /msu7/msu7.genome_tran
```
&emsp; 2. Mapping by HISAT2

```
  # Check sequencing quality 
  find -L /RNAseq/00_fastq -name "*.fq.gz" | \
      xargs -n 1 -P 4 -I PREFIX \
      sh -c '
          sample=`basename PREFIX | cut -d "." -f 1`

          echo "[`date`]: Start mapping ${sample} ... "
          fastqc --noextract -o /RNAseq/00_fastq/stats --format fastq \
              --threads 4 PREFIX
      '

  # Align reads and sort bam files
  find -L /RNAseq/00_fastq/ -name "*_1.clean.fq.gz" | sed 's/_1.clean.fq.gz$//' | \
      xargs -n 1 -P 5 -I PREFIX \
      sh -c '

          sample=`basename PREFIX | cut -d "." -f 1`

          echo "[`date`]: Start mapping ${sample} ... "

          read1=PREFIX"_1.clean.fq.gz"
          read2=PREFIX"_2.clean.fq.gz"

          ## Align reads with hisat2
          hisat2 -p 4 --dta -x /msu7/msu7.genome_tran  \
              -1 ${read1} -2 ${read2} -S /RNAseq/01_assembly/${sample}.msu7.hisat2.sam \
              > /RNAseq/01_assembly/${sample}.msu7.hisat2.log 2>&1

          ## sort bam file
          java -Djava.io.tmpdir=/tmp -jar picard-2.9.0/picard.jar SortSam \
              I=/RNAseq/01_assembly/${sample}.msu7.hisat2.sam \
              O=/RNAseq/01_assembly/${sample}.msu7.hisat2.sort.bam \
              SORT_ORDER=coordinate \
              >> /RNAseq/01_assembly/${sample}.msu7.hisat2.log 2>&1 && \
              rm -v /RNAseq/01_assembly/${sample}.msu7.hisat2.sam
          ## index bam file
          samtools index /RNAseq/01_assembly/${sample}.msu7.hisat2.sort.bam
      '
```
* Estimate abundance for each sample by StringTie

```
  # Count abundance and create table counts
  mkdir -p /RNAseq/04_count
  find /RNAseq/01_assembly/ -name "*msu7.hisat2.sort.bam" | \
      xargs -n 1 -P 4 -I PREFIX \
      sh -c '
          sample=`basename PREFIX | cut -d "." -f 1`

          echo ${sample}

          mkdir -p /RNAseq/04_count/${sample}/
          stringtie-1.3.4b.Linux_x86_64/stringtie \
              -e -B -p 4 \
              -G /MSU_osa1r7/all.mod.gtf \
              -o /RNAseq/04_count/${sample}/${sample}.msu7.hisat2.sort.trans.gtf \
              /RNAseq/01_assembly/${sample}.msu7.hisat2.sort.bam
      '
  # Extract FPKM, TPM and coverage for each gene
  find /RNA_seq/04_count/  -name "*.sort.trans.gtf" | \
      xargs -n 1 -P 6 -I PREFIX \
      sh -c '
          sample=`basename PREFIX | cut -d "." -f 1`
          echo "[`date`]: Start processing ${sample} ... "
          perl count2table.pl \
              -i PREFIX \
              > PREFIX.tab
      '
```
* Identify DEGs by ballgown  
```
## assign experimental design
perl experimental_design.pl \
    -i experimental_design.txt \
    -cd /RNA_seq/04_count \
    -od RNA_seq/05_deg
# Phenotype name for control sample in "experimental_design.txt" (tab-separated) must be first order in alphabetical order.

## ballgown installation
# source("http://bioconductor.org/biocLite.R")
# biocLite("ballgown")
# biocLite("RSkittleBrewer")
# biocLite("devtools")
# devtools::install_github('alyssafrazee/RSkittleBrewer', force = TRUE)

```


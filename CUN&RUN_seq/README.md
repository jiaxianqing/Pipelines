## CUT&RUN-seq analysis

> Dr. Xianqing Jia (jiaxianqing@nwu.edu.cn)   
> Latest update: 2024-01-29   
> Lab website: https://jialab.life/   

This pipeline is based on:
- deepTools: https://deeptools.readthedocs.io/en/develop/index.html
- MACS3: https://github.com/macs3-project/MACS
- ChIPseeker: https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html

It includes:


---
### 1. Trim fastq files with fastp
```
    mkdir -p /mnt/hdd1/Data/rice/projects/WRKY45/00_fastq/clean_trimed
    find -L /mnt/hdd1/Data/rice/projects/WRKY45/00_fastq/ -name "*_1.clean.fq.gz" | \
        sed 's/_1.clean.fq.gz$//' | xargs -n 1 -P 2 -I PREFIX \
        sh -c '
            
            sample=`basename PREFIX | cut -d "_" -f 1`
            
            echo "[`date`]: Start processing ${sample} ... "
            
            read1=PREFIX"_1.clean.fq.gz"
            read2=PREFIX"_2.clean.fq.gz"

            /mnt/hdd1/Data/biosoft/Fastp/fastp \
                -i ${read1} -I ${read2} \
                -o /mnt/hdd1/Data/rice/projects/WRKY45/00_fastq/clean_trimed/${sample}_1.trimed.fq.gz \
                -O /mnt/hdd1/Data/rice/projects/WRKY45/00_fastq/clean_trimed/${sample}_2.trimed.fq.gz \
                --html /mnt/hdd1/Data/rice/projects/WRKY45/00_fastq/clean_trimed/${sample}.trimed.html \
                --json /mnt/hdd1/Data/rice/projects/WRKY45/00_fastq/clean_trimed/${sample}.trimed.json \
                > /mnt/hdd1/Data/rice/projects/WRKY45/00_fastq/clean_trimed/${sample}.trimed.log 2>&1
        '

```
### 2. Mapping

```
mkdir -p /mnt/hdd1/Data/rice/projects/WRKY45/01_assembly2
find -L /mnt/hdd1/Data/rice/projects/WRKY45/00_fastq/clean_trimed/ -name "*_1.trimed.fq.gz" | sed 's/_1.trimed.fq.gz$//' | \
    xargs -n 1 -P 2 -I PREFIX \
    sh -c '
        
        sample=`basename PREFIX | cut -d "_" -f 1`
        
        echo "[`date`]: Start mapping ${sample} ... "
        
        read1=PREFIX"_1.trimed.fq.gz"
        read2=PREFIX"_2.trimed.fq.gz"
        
        ## Align reads with bowtie2, --very-sensitive
        /mnt/hdd1/Data/biosoft/bowtie2-2.4.4-linux-x86_64/bowtie2 -p 3 --very-sensitive \
            -1 ${read1} -2 ${read2} \
            -x /mnt/hdd1/Data/rice/ref/bowtie2_index/IRGSP \
            -S /mnt/hdd1/Data/rice/projects/WRKY45/01_assembly2/${sample}.IRGSP.bowtie2.sam \
            > /mnt/hdd1/Data/rice/projects/WRKY45/01_assembly2/${sample}.IRGSP.bowtie2.log 2>&1
        
        ## sort bam file
        java -Djava.io.tmpdir=/mnt/hdd1/Data/tmp -jar /mnt/hdd1/Data/biosoft/picard/picard-2.23.3.jar SortSam \
            I=/mnt/hdd1/Data/rice/projects/WRKY45/01_assembly2/${sample}.IRGSP.bowtie2.sam \
            O=/mnt/hdd1/Data/rice/projects/WRKY45/01_assembly2/${sample}.IRGSP.bowtie2.sort.bam \
            SORT_ORDER=coordinate \
            >> /mnt/hdd1/Data/rice/projects/WRKY45/01_assembly2/${sample}.IRGSP.bowtie2.log 2>&1 && \
            rm -v /mnt/hdd1/Data/rice/projects/WRKY45/01_assembly2/${sample}.IRGSP.bowtie2.sam
        ## index bam file
        samtools index /mnt/hdd1/Data/rice/projects/WRKY45/01_assembly2/${sample}.IRGSP.bowtie2.sort.bam
    '
```


### 3. filter reads
(MQ<10, Length<140, remove duplicates)
```
find /mnt/hdd1/Data/rice/projects/WRKY45/01_assembly2/ -name "*.bam" | \
    xargs -n 1 -P 2 -I PREFIX \
    sh -c '
        
        sample=`basename PREFIX | cut -d "." -f 1`
        
        echo "[`date`]: Start dedup ${sample} ... "
        ## filt bam file --ignoreDuplicates
        alignmentSieve --minMappingQuality 10 --minFragmentLength 140 --ignoreDuplicates --numberOfProcessors 2 \
            -b /mnt/hdd1/Data/rice/projects/WRKY45/01_assembly2/${sample}.IRGSP.bowtie2.sort.bam \
            -o /mnt/hdd1/Data/rice/projects/WRKY45/01_assembly2/${sample}.IRGSP.bowtie2.sort.filtered.bam \
            --filterMetrics /mnt/hdd1/Data/rice/projects/WRKY45/01_assembly/45-flag-9-1.IRGSP.bowtie2.sort.filtered.log
        
        ## sort bam file
        java -Djava.io.tmpdir=/mnt/hdd1/Data/tmp -jar /mnt/hdd1/Data/biosoft/picard/picard-2.23.3.jar SortSam \
            I=/mnt/hdd1/Data/rice/projects/WRKY45/01_assembly2/${sample}.IRGSP.bowtie2.sort.filtered.bam \
            O=/mnt/hdd1/Data/rice/projects/WRKY45/01_assembly2/${sample}.IRGSP.bowtie2.sort.filtered.sort.bam \
            SORT_ORDER=coordinate \
            >> /mnt/hdd1/Data/rice/projects/WRKY45/01_assembly2/${sample}.IRGSP.bowtie2.log 2>&1 && \
            rm -v /mnt/hdd1/Data/rice/projects/WRKY45/01_assembly2/${sample}.IRGSP.bowtie2.sort.filtered.bam
        
        ## index bam file
        samtools index /mnt/hdd1/Data/rice/projects/WRKY45/01_assembly2/${sample}.IRGSP.bowtie2.sort.filtered.sort.bam
    '
```

### 3. peak calling

#### 1) Effective Genome Size
https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html

```
    # pip install khmer
    unique-kmers.py -k 150 /mnt/hdd1/Data/rice/ref/IRGSP-1.0_genome.fasta \
        > /mnt/hdd1/Data/rice/ref/IRGSP-1.0_genome.fasta.unique-kmers.txt
    # rice genome: 355020671
```

#### 2) calling peaks
```
    mkdir -p /mnt/hdd1/Data/rice/projects/WRKY45/02_peak2/

    python3 -m venv MyPythonEnv/
    source MyPythonEnv/bin/activate

    ## per sample
    # "-g" is mappable genome size or effective genome size
    for file in `ls /mnt/hdd1/Data/rice/projects/WRKY45/01_assembly2/*.IRGSP.bowtie2.sort.filtered.sort.bam`;
    do
        samp=`basename ${file} | cut -f 1 -d "."`
        echo ${samp}
        macs3 callpeak -f BAM -g 355020671 -B -q 0.05 \
            --nomodel --extsize 200 --shift -100 --keep-dup all --call-summits -n ${samp} \
            -t ${file} \
            --outdir /mnt/hdd1/Data/rice/projects/WRKY45/02_peak2/
    done

    ## control vs. wrky45
    macs3 callpeak -f BAM -g 374471240 -B -q 0.05 \
        --nomodel --extsize 200 --shift -100 --keep-dup all --call-summits -n WRKY45 \
        -c /mnt/hdd1/Data/rice/projects/WRKY45/01_assembly2/45-flag-9-ck.IRGSP.bowtie2.sort.filtered.sort.bam \
        -t /mnt/hdd1/Data/rice/projects/WRKY45/01_assembly2/45-flag-9-1.IRGSP.bowtie2.sort.filtered.sort.bam /mnt/hdd1/Data/rice/projects/WRKY45/01_assembly2/45-flag-9-2.IRGSP.bowtie2.sort.filtered.sort.bam \
        --outdir /mnt/hdd1/Data/rice/projects/WRKY45/02_peak2/
```

### 4. annotation

```
BiocManager::install("ChIPseeker")
BiocManager::install("GenomicFeatures")
```

```
library(ChIPseeker)
library(GenomicFeatures)
spompe <- makeTxDbFromGFF('all.mod.gtf')

# 03
peak <- readPeakFile('WRKY76_peaks.xls')
peakAnno <- annotatePeak(peak, tssRegion = c(-3000, 3000), TxDb = spompe)
write.table(peakAnno, file = 'WRKY76_peaks.info.xls',sep = '\t', quote = FALSE, row.names = FALSE)

# 04
plotAnnoBar(peakAnno)
vennpie(peakAnno)
plotAnnoPie(peakAnno)

# 05
peak1 <- readPeakFile('Xu_MUT_rep1_summits.bed')
peak2 <- readPeakFile('Xu_MUT_rep2_summits.bed')
peaks <- list(peak1 = peak1, peak2 = peak2)
peakAnnoList <- lapply(peaks, annotatePeak, TxDb = spompe, tssRegion = c(-3000, 3000), addFlankGeneInfo = TRUE, flankDistance = 5000)
write.table(peakAnnoList[1], file = 'peak1.txt',sep = '\t', quote = FALSE, row.names = FALSE)
write.table(peakAnnoList[2], file = 'peak2.txt',sep = '\t', quote = FALSE, row.names = FALSE)

# 06.
plotAnnoBar(peakAnnoList)
vennpie(peakAnnoList[[1]])
plotDistToTSS(peakAnnoList)
```

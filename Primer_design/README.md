## Pipeline for batch primer design
If you want to design many primer pairs for different sites along the genome, this pipeline might be a good choice.


## Software dependencies
Primer3 <https://github.com/primer3-org/primer3><br>
BLAST <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST>

## Steps

#### Step1. Prepare sequences and input files

* Extract seqeunces from the whole genome sequence file. Please skip this, if you already have a sequence file in fasta style. bedtools is used to extact sequences. An example is given here:
```
awk 'BEGIN{OFS="\t"} {print $1,$2-1000,$2+1000,$3;}' candidate.site.list | \
    bedtools getfasta -name -bed - -fi IRGSP_1.fa -fo candidate.site.list.ex1k.fa
```
Note: **candidate.site.list** file has four rows: chromosome,  start, end, name, an example:  
>chr01   1000    2000    test1<br>
>chr03   3000    4000    test2<br>

* Transform FASTA file to input file for Primer3.
```
perl fasta2primer3.pl --task pick_pcr_primers --target_region 900,200 --size-region 400-700 --number 30 \
    --input candidate.site.list.ex1k.fa \
    > candidate.site.list.ex1k.fa.primer3.input
```
#### Step2. Run Primer3

```
/home/jxq/biosoft/primer3-2.4.0/src/primer3_core -format_output \
    -p3_settings_file=primer3web_v4_0_0_default_settings.mod.txt \
    -error=candidate.site.list.ex1k.fa.primer3.input.err \
    -output=candidate.site.list.ex1k.fa.primer3.output \
    candidate.site.list.ex1k.fa.primer3.input
```

#### Step3. Find optimal primers

```
perl /home/jxq/Data/scripts/my_scripts/primer3_to_blast.pl \
    -i /home/jxq/Data/rice/project/head_date2/07_genes/primer/candidate.site.list.ex1k.fa.primer3.output \
    -s /home/jxq/Data/rice/project/head_date2/07_genes/primer/candidate.site.list.ex1k.fa.primer3.output.stats \
    -o /home/jxq/Data/rice/project/head_date2/07_genes/primer/candidate.site.list.ex1k.fa.primer3.output.fas
```



```
mkdir -p /home/jxq/Data/rice/blast/library
makeblastdb -dbtype nucl -in /home/jxq/Data/rice/ref/IRGSP-1.0_genome.fasta \
    -out /home/jxq/Data/rice/blast/library/IRGSP.genome \
    -logfile /home/jxq/Data/rice/blast/library/IRGSP.genome.log
```

```
blastn -query /home/jxq/Data/rice/project/head_date2/07_genes/primer/candidate.site.list.ex1k.fa.primer3.output.fas \
    -db /home/jxq/Data/rice/blast/library/IRGSP.genome \
    -task blastn -evalue 1 -num_threads 4 -max_target_seqs 10 -max_hsps 10 \
    -outfmt '7 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend qcovs sstart send evalue bitscore' \
    -out /home/jxq/Data/rice/project/head_date2/07_genes/primer/candidate.site.list.ex1k.fa.primer3.output.fas.IRGSP.genome.blastn
```

```
perl /home/jxq/Data/scripts/my_scripts/blast_stats.pl \
    -i /home/jxq/Data/arabidopsis/projects/seed_traits/04_primer/cds_last300.bed.fa.primer3.output.fas.TAIR10_genome.blastn \
    -o /home/jxq/Data/arabidopsis/projects/seed_traits/04_primer/cds_last300.bed.fa.primer3.output.fas.TAIR10_genome.blastn.stats
```

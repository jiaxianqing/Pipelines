## Pipeline for batch primer design
If you want to design many primer pairs for different sites along the genome, this pipeline might be a good choice.


## Software dependencies
Primer3 <https://github.com/primer3-org/primer3><br>
BLAST <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>

## Steps

#### Step1. Prepare sequences and input files

* Extract seqeunces from the whole genome sequence file. Please skip this, if you already have a sequence file in fasta style. bedtools is used to extact sequences. An example is given here:

        awk 'BEGIN{OFS="\t"} {print $1,$2-1000,$2+1000,$3;}' candidate.site.list | \
            bedtools getfasta -name -bed - -fi IRGSP_1.fa -fo candidate.site.list.ex1k.fa

Note: **candidate.site.list** file has four rows: chromosome,  start, end, name, an example:
>chr01   1000    2000    test1<br>
>chr03   3000    4000    test2<br>


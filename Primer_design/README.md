## Pipeline for batch primer design
If you want to design many primer pairs for different sites along the genome, this pipeline might be a good choice.


## Software dependencies
Primer3 <https://github.com/primer3-org/primer3>

BLAST <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>

## Steps

#### Step1. Prepare sequences and input files

* Extract seqeunces from the whole genome sequence file. Please skip this, if you already have a sequence file in fasta style. bedtools is used to extact sequences. Give a example here:
    awk 'BEGIN{OFS="\t"} {print $1,$2-1000,$2+1000,$3;}' candidate.site.list | \
        getfasta -name -bed - -fo candidate.site.list.ex1k.fa


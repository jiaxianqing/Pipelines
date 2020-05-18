## SRA process   
This pipeline is for SRA processing. Including SRA download, SRA2FASTQ and how to upload raw data.

---

* **SRA download**  
We could download SRA files from NCBI database by FTP links (download quickly). An example, for SRR3056114, ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR305/SRR3056114/SRR3056114.sra. Here is a more general way to download SRA files, by [SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit).

1. On Windows  
The script is wrapped in a bat file `sra_download.bat`. One can directly run `sra_download.bat` and enter the needed informations as prompted.  
2. On Linux  
After installation of SRA-Toolkit for Linux, one can use `prefetch` to download files.

* **SRA files to FASTQ files**  

```
# file check
./sratoolkit.2.9.0-centos_linux64/bin/vdb-validate SRR3056114.sra

# sra2fastq
/sratoolkit.2.9.0-centos_linux64/bin/fastq-dump --split-files --gzip \
  --outdir ./fastq SRR3056114.sra
```

* **Upload raw data to NCBI SRA database**  


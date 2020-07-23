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
Usually, we need to upload our raw data to some public database, like NCBI. But for NCBI, there is no a detailed description or protocol for "how to upload raw data". Here an example for upload DNA sequencing data to NCBI is given:

1. Create a new `BioProject`  
Assume that we already have an account (if not, please [create one](https://submit.ncbi.nlm.nih.gov/subs/)). Then please go to [NCBI Sign In Page](https://submit.ncbi.nlm.nih.gov/subs/) and login. Next step is to create a new `BioProject`.


<img src="https://github.com/jiaxianqing/Pipelines/blob/master/SRA_process/new_bioproject.png" width = "30%" height = "30%" div align = "center" />

2. Create a new `BioSample`  
Visit https://submit.ncbi.nlm.nih.gov/subs/biosample/ and create a new BioSample.  

<img src="https://github.com/jiaxianqing/Pipelines/blob/master/SRA_process/new_biosample.png" width = "30%" height = "30%" div align = "center" />

Then we can download template table for BioSample info.  

<img src="https://github.com/jiaxianqing/Pipelines/blob/master/SRA_process/new_biosample2.png" width = "50%" height = "50%" div align = "center" />

Here is an example for BioSample info:

<img src="https://github.com/jiaxianqing/Pipelines/blob/master/SRA_process/Sanple_info.png" width = "50%" height = "50%" div align = "center" />

Then we need to go back to submission list to download BioSample IDs.

<img src="https://github.com/jiaxianqing/Pipelines/blob/master/SRA_process/biosample_accessions.png" width = "50%" height = "50%" div align = "center" />

3. Creat a new `preload`  
This step is to upload your files to database. Go to https://submit.ncbi.nlm.nih.gov/subs/sra/.

<img src="https://github.com/jiaxianqing/Pipelines/blob/master/SRA_process/new_preload.png" width = "50%" height = "50%" div align = "center" />

Please choose `FTP upload` and we will got Username and Password.

<img src="https://github.com/jiaxianqing/Pipelines/blob/master/SRA_process/file_upload2.jpg" width = "50%" height = "50%" div align = "center" />

Now, it's time to upload our files! I recommend to use the FTP tool, [FileZilla Client](https://filezilla-project.org/)  

<img src="https://github.com/jiaxianqing/Pipelines/blob/master/SRA_process/filezilla.png" width = "50%" height = "50%" div align = "center" />  
<img src="https://github.com/jiaxianqing/Pipelines/blob/master/SRA_process/filezilla3.jpg" width = "50%" height = "50%" div align = "center" />

4. Create a new `SRA submission`  
Go to https://submit.ncbi.nlm.nih.gov/subs/sra/ and select `New submission` to upload `SRA_metadata` information.  


<img src="https://github.com/jiaxianqing/Pipelines/blob/master/SRA_process/new_sra_submission.png" width = "30%" height = "30%" div align = "center" />  
<img src="https://github.com/jiaxianqing/Pipelines/blob/master/SRA_process/sra_submission_info.png" width = "50%" height = "50%" div align = "center" />

---
Create: 2020-05-18  
Update: 2020-05-19  
Xianqing Jia  
jiaxq.nju@gmail.com  

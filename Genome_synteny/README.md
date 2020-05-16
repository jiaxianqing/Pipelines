## Genome synteny analysis

This pipeline is followed [MCscan (Python version)](https://github.com/tanghaibao/jcvi/wiki/MCscan-(Python-version)).

---

### Software dependencies  
[LAST](http://last.cbrc.jp/) and [jcvi](https://github.com/tanghaibao/jcvi/releases)

### Pipeline  
* **Installation**  

&emsp; 1. LAST
```
  unzip last-1060.zip
  cd last-1060
  make
  
  # Add scripts and src directories to $PATH
```
&emsp; 2. jcvi
```
  pip install jcvi

  # Possible Error reporting and resolution:
  # E1: AttributeError: module 'pandas' has no attribute 'core'
  # R1: /Python-3.6.1/bin/pip3 install dask --upgrade
  #     /Python-3.6.1/bin/pip3 install pandas==0.22 --user
```
* **Data preparetion**  

CDS sequences file (`FASTA` format) and coordinates file (`BED` format) are needed for MCscan. An example for `BED` file (tab-separated):
>Chr1	2902	10817	LOC_Os01g01010.1	0	+  
>Chr1	2983	10562	LOC_Os01g01010.2	0	+  
>Chr1	11217	12435	LOC_Os01g01019.1	0	+  
>Chr1	12647	15915	LOC_Os01g01030.1	0	+  

For several frequently-used plant genomes, we can conveniently download these needed data from [Phytozome](https://phytozome.jgi.doe.gov/pz/portal.html) by jcvi.
```
  # Login phytozome account
  $ python3 -m jcvi.apps.fetch phytozome
  Phytozome Login: xxxxxxxx
  Phytozome Password:

  # Download data from Phytozome
  python3 -m jcvi.apps.fetch phytozome Osativa,Atrichopoda,Spolyrhiza,Vvinifera

  # GFF to BED and convert CDS file
  for file in `ls ./*.gene.gff3.gz`
  do
      echo ${file}
      name=`basename ${file} | cut -f 1 -d "_"`
      python3 -m jcvi.formats.gff bed --type=mRNA --key=Name \
          ${file} \
          -o ${name}.bed
      python3 -m jcvi.formats.fasta format \
          ${file} \
          ${name}.cds
  done
```
* **Run pairwise synteny search**

```
  # For example, compare Osativa to others
  for file in `ls ./*.bed`
  do
      name=`basename ${file} | cut -f 1 -d "."`
      if [ ${name} == "Osativa" ]
      then
          continue
      fi
      echo ${name}
      python3 -m jcvi.compara.catalog ortholog Osativa ${name}
  done
```
* **Synteny visualization**  

After pairwise synteny search, there are two strategies to show resluts, macrosynteny (genomic) and microsynteny (certain region). One can find minces work flow at https://github.com/tanghaibao/jcvi/wiki/MCscan-(Python-version). Here a example for microsynteny is given:  
```
# Construct multi-synteny blocks
for file in `ls ./*.lifted.anchors`
do
    name=`basename ${file} | cut -f 2 -d "."`
    echo ${name}
    python3 -m jcvi.compara.synteny mcscan \
        Osativa.bed \
        ${file} \
        --iter=5 \
        -o Osativa.${name}.i5.blocks
done

# Merge synteny blocks
python3 -m jcvi.formats.base join \
    Osativa.Atrichopoda.i1.blocks \
    Osativa.Spolyrhiza.i1.blocks \
    Osativa.Vvinifera.i1.blocks --noheader | \
        cut -f 1,2,4,6,8 > Osativa2all.i1.blocks

# Merge bed files
cat Osativa.bed Atrichopoda.bed Spolyrhiza.bed Vvinifera.bed \
    > 4species.bed

# Plot synteny image
python3 -m jcvi.graphics.synteny \
    test.blocks \
    4species.bed \
    test.layout

# "test.blocks" is captured from "Osativa2all.i1.blocks"
# "test.layout" file (blocks.layout file) (comma+space-separated):
  #x, y, rotation, ha, va, color, ratio, label
  0.5, 0.3, 0, left, center, , 1, Osativa
  0.5, 0.4, 0, left, center, , 1, Atrichopoda
  0.5, 0.5, 0, left, center, , 1, Spolyrhiza
  0.5, 0.6, 0, left, center, , 1, Vvinifera
  # edges
  e, 0, 1
  e, 1, 2
  e, 2, 3
```
The resulting plot below:

<img src="https://github.com/jiaxianqing/Pipelines/blob/master/Genome_synteny/genome_synteny.png"  div align = "center" width="75%" height="75%" />


---

* **Plot synteny image**  
Assume that you already have synteny list of genes, then if you want to show their synteny relationship.

```
  cat Chr2.bed Chr6.bed \
      > Chr2_Chr6.bed

  /home/jxq/bin/Python-3.6.1/bin/python3 -m jcvi.graphics.synteny \
      Chr2_Chr6.blocks \
      Chr2_Chr6.bed \
      Chr2_Chr6.layout
```

<img src="https://github.com/jiaxianqing/Pipelines/blob/master/Genome_synteny/Chr2_Chr6.png"  div align = "center" width="75%" height="75%" />

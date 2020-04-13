## Summary
Here three common strategies used in phylogeny construction are introduced.
* [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/index.html)
* [FastTreeMP](http://www.microbesonline.org/fasttree/)
* [IQ-TREE](http://www.iqtree.org/)

### Step 1. Sequence alignment
Three widely-used softwares and their usages are listed below.
* [MAFFT](https://mafft.cbrc.jp/alignment/software/)
```
  mafft --thread 8 --maxiterate 1000 --localpair os-cr.fas \
      > os-cr.mafft.aln
```
* [Clustalo](http://www.clustal.org/omega/)
```
  clustalo --threads=8 -v -i os-cr.fas \
      -o os-cr.clustalo.aln
```
* [MUSCLE](http://www.drive5.com/software.html)
```
  muscle -in os-cr.fas -out os-cr.muscle.aln -quiet
```
Note: "Phylip" input format is required for some softwares, like [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/index.html) and [PhyML](http://www.atgc-montpellier.fr/phyml/). Format conversion is present by [fasta2relaxedPhylip.pl](https://github.com/npchar/Phylogenomic)
```
  perl fasta2relaxedPhylip.pl -f os-cr.mafft.aln \
      -o os-cr.mafft.aln.phy
```

## Step 2. Phylogeny construction
* RAxML
```
raxmlHPC-PTHREADS-SSE3 \
    -s os-cr.mafft.aln.phy \
    -m PROTGAMMAVTX -f a -x 12345 -p 12345 -# 1000 -n tre -T 12 \
    -w ./tree
```
* FastTree (recommended)
```
  FastTreeMP \
      -log os-cr.mafft.aln.FastTree.log \
      os-cr.mafft.aln \
      > os-cr.mafft.aln.FastTree.nwk
```
* IQ-TREE (recommended)
```
  iqtree -s os-cr.mafft.aln \
      -m MFP -nt AUTO -ntmax 10 -bb 1000 -alrt 1000 \
      -pre os-cr.mafft.aln.iqtree
```

### Step 3. Tree editing 
Final, we could edit or anotate phylogenetic tree by [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) or [iTOL](https://itol.embl.de/index.shtml) (recommended). Here is a final version of my "Dof tree".

<img src="https://github.com/jiaxianqing/Pipelines/raw/master/Phylogeny_construction/Dof%20genes.jpg" width = "500" height = "500" div align = "center" />

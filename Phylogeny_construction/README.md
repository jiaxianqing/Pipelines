## Phylogeny construction


---

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


## Species phylogeny construction
Species phylogeny is usually constructed by genome-wide concatenated orthologous genes (single-copy genes or low-copy genes (one to two copies for each genome)). [OrthoFinder](https://github.com/davidemms/OrthoFinder) is used to identify gene orthologs. Multiple sequence alignment was performed by MAFFT. The tree was built by RAxML or IQ-TREE.

---

* **Identify gene orthologs**
```
  # Assume that you already put whole-genome protein file for each genome into ./fasta directory
  orthofinder -t 12 -a 8 -og \
      -S diamond \
      -p ./tmp \
      -f ./fasta \
      > ./fasta/orthofinder.diamond.process.log 2>&1
# Single copy orthogroups are listed in "SingleCopyOrthogroups.txt".

# Find low copy genes and select one of them in each orthogroups each species.
perl find_low_copy_genes.pl -- genecount Orthogroups.GeneCount.csv \
    --orthogroups Orthogroups.csv \
    --lowcopylist Orthogroups.GeneCount.low_copy.list \
    --lowcopygroup Orthogroups.GeneCount.low_copy.gene_list.csv
```
* **Multiple sequence alignment**  
There are two strategies for species phylogeny construction by single- or low- copy genes, `concatenation` and `coalescence`.   
1. `Concatenation`: each single copy gene of different species was isolated for multi-sequence alignment, and then the single copy genes were concatenated end to end into a supergene matrix, which was finally used to construct a phylogenetic tree.  
2. `Coalescence`: multiple sequence alignment was carried out for each single copy gene among different species, and the gene tree corresponding to each single copy gene was constructed. Then, the gene trees corresponding to all single copy genes were combined to reconstruct the corresponding species tree. 
Here, this pipeline is based on `Concatenation` way:
```
  # Extract sequences from each speice and put them into ./Orthogroups/
  mkdir ./Orthogroups
  perl find_gene_seq.pl \
      --seqdir ./fasta --suffix fa \
      --lowcopygroup Orthogroups.GeneCount.low_copy.gene_list.csv \
      --outdir ./Orthogroups
  #If some sequences are lost, please try: dos2unix [workdir]
  
  # Multiple sequence alignment
  for file in `ls ./Orthogroups/*.fa`
  do
      echo ${file}
      mafft --thread 8 --maxiterate 1000 --localpair  ${file} > ${file}.aln
  done
  
  # Merge aligned sequences (Concatenation)
  perl merge_aln_seq.pl \
      --indir ./Orthogroups --suffix aln \
      --lowcopygroup Orthogroups.GeneCount.low_copy.gene_list.csv \
      --output low_copy_gene.fa.maerge.aln
```
* **Construct species phylogeny**
```
# RAxML
perl fasta2relaxedPhylip.pl -f low_copy_gene.fa.maerge.aln \
    -o low_copy_gene.fa.maerge.aln.phy
raxmlHPC-PTHREADS-SSE3 \
    -s low_copy_gene.fa.maerge.aln.phy \
    -m PROTGAMMAJTTX -f a -x 12345 -p 12345 -N 1000 -n tre -T 12 \
    -w ./

# FastTree
FastTreeMP \
    -log low_copy_gene.fa.maerge.aln.FastTree.log \
    low_copy_gene.fa.maerge.aln \
    > low_copy_gene.fa.maerge.aln.FastTree.nwk

# IQ-TREE (recommended)
iqtree \
    -s low_copy_gene.fa.maerge.aln \
    -m MFP -nt AUTO -ntmax 16 -bb 1000 -alrt 1000 -asr \
    -pre low_copy_gene.fa.maerge.aln.iqtree
```

An example for final visualization:  

<img src="https://github.com/jiaxianqing/Pipelines/blob/master/Phylogeny_construction/species_tree.png" width = "25%" height = "25%" div align = "center" />

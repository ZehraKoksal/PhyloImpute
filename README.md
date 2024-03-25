# PhyloImpute
Impute missing data on non-recombining DNA using SNP phylogeny

### 1) About PhyloImpute
PhyloImpute complements SNP data on non-recombining DNA, such as the human Y chromosome, by leveraging the SNPs' phylogeny as reported in a phylogenetic tree.

### 2) Installation
Operating system: Linux

Type in the shell:
```
git clone https://github.com/ZehraKoksal/PhyloImpute.git
cd PhyloImpute/
python PhyloImpute.py -h
```

### 3) Algorithm and Commands
PhyloImpute imputes missing data by assuming that the SNPs in a clade of the phylogenetic tree leading up to a SNP with a derived allele are derived as well. SNPs on parallel branches are expected to be ancestral.

PhyloImpute requires an input file (more details below), the desired output file and the phylogenetic tree. 

#### 3.1) Input file
The user is required to provide the path to the input file in the tab-separated _.csv_ format. 

<img src="/test_run/images/Input.png" alt="Input file style" width="700"/>

The rows present variants, the columns individuals.
The header row should present the individuals' labels (blue) and the second column the variant names (orange). The table contains the oberseved allelic states (green) which can be ancestral **A**, derived **D**  or missing **X** for each variant.



#### 3.2) Phylogenetic tree
The phylogenetic tree needs to contain at least one of the SNPs in the input file in the exact same nomenclature. 

Currently, a pre-processed phylogenetic tree is available for the human Y chromosome (Minimal Y tree):
```
python PhyloImpute.py input_file.csv/ output_file/ --tree Y_minimal
```

Additionally, custom phylogenetic trees can be provided:
```
python PhyloImpute.py input_file.csv/ output_file/ --customtree example_custom_tree_Minimal_Y_tree
```

The custom phylogenetic trees need to be made available by the user in the tab-separated _.csv_ format. One example can be viewed in the test_run folder provided here. It follows the structure presented here and matching the ISOGG tree nomenclature (https://isogg.org/tree/):

<img src="/test_run/images/Custom_tree.png" alt="Input file style" width="700"/>


#### 3.3) Output file
The outcome are the observed (**D**,**A**,**X**) and imputed (**d**,**a**,**X**) allelic states for the initially reported SNPs complemented with the SNPs in the phylogenetic tree. 

<img src="/test_run/images/Output_partly.png" alt="Input file style" width="700"/>

PhyloImpute cross-references the allelic states of all observed SNPs with the SNP relationships in the phylogenetic tree to verify the accuracy of the phylogenetic tree and the sequencing data. SNPs that cause contradictions between the provided input data and the imputations from the phylogenetic tree, keep the observed alleles and are stored in an additional file **questionable_SNPs.csv** for manual inspection.
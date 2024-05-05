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

PhyloImpute requires an input file, the desired output file, and a phylogenetic tree. 

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

The custom phylogenetic trees need to be made available by the user in the tab-separated _.csv_ format. One example can be viewed in the test_run folder provided here. It follows the structure presented here and matches the ISOGG tree nomenclature (https://isogg.org/tree/):

<img src="/test_run/images/Custom_tree_hg.png" alt="Input file style" width="700"/>

SNPs that cannot be separated ("equal") are divided by commas in the same branch (green). SNPs of downstream branches are presented in the row below with one additional indentation using a tab (orange). And SNPs from parallel branches are on separate, mutually exclusive branches (blue).

#### 3.3) Output file
The outcome are the observed (**D**,**A**,**X**) and imputed (**d**,**a**,**X**) allelic states for the initially reported SNPs complemented with the SNPs in the phylogenetic tree and a rooting SNP (ROOT). 

<img src="/test_run/images/Output_partly.png" alt="Input file style" width="400"/>

PhyloImpute cross-references the allelic states of all observed SNPs with the SNP relationships in the phylogenetic tree to verify the accuracy of the phylogenetic tree and the sequencing data. 

Following the tree branch with the most SNPs that are present in the analyzed sequence ("maximum resolution tree branch"), the haplogroup of the sample will be predicted. The output is given in the file **haplogroups.csv**, where each row present a sample, the predicted haplogroup and two statistical values to assess the reliability of the prediction. The first value is the **confidence value**, which counts the different derived alleles of SNPs from the **maximum resolution tree** branch divided by the total variants in that branch. A low **confidence value** can be simply due to low sequence coverage and does not necessarily mean that the predicted haplogroup is wrong. The second value, the **penalty value**, counts the different ancestral alleles of variants from the **maximum resolution tree** branch divided by the total variants in that branch. While ancestral alleles could be the consequence of backmutations, a high penalty value generally hints towards an incorrect haplogroup prediction. 

Contradictions can be observed for SNPs that have different allelic states between the provided input data and the imputations from the phylogenetic tree. In these cases, PhyloImpute keeps the observed allelic states and stores the contradicting SNPs in an additional file **conflicting_SNPs.csv** for manual inspection of the phylogenetic tree and sequencing data.
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

PhyloImpute can be run with a pre-processed csv input file or vcf files:

#
#### 3.1) CSV file
```
python PhyloImpute.py -input_format csv -input ./test_run/testdata.csv -output ./output -tree Y_minimal

```
**Parameters:**

**-input_format** csv

**-input** path to the input file

**-tree** path to the available phylogenetic tree {Y_minimal} [is mutually exclusive with -customtree] 

**-customtree** path to custom phylogenetic tree [is mutually exclusive with -tree] 

**-output** path to the existing folder for the output files



#### 3.1.1) CSV Input file
The user is required to provide the path to the input file in the tab-separated _.csv_ format. 

<img src="/test_run/images/input_csv.png" alt="Input file style" width="350"/>

The rows present variants, the columns individuals.
The header row should present the individuals' labels (blue) and the second column the variant names (orange). The table contains the oberseved allelic states (green) which can be ancestral **A**, derived **D**  or missing **X** for each variant.

#### 3.1.2) Phylogenetic tree
The phylogenetic tree needs to contain at least one of the SNPs in the input file in the exact same nomenclature. 

#### 3.1.2.1) Pre-processed phylogenetic tree
Currently, a pre-processed phylogenetic tree is available for the human Y chromosome (Minimal Y tree):
```
python PhyloImpute.py -input_format csv -input ./test_run/testdata.csv -output ./output -tree Y_minimal
```

#### 3.1.2.2) Custom phylogenetic tree
Alternatively, custom phylogenetic trees can be provided:
```
python PhyloImpute.py -input_format csv -input ./test_run/testdata.csv -output ./output -customtree ./test_run/minimal_y_tree_hgs_custom_example.csv
```

The custom phylogenetic tree need to be made available by the user in the tab-separated _.csv_ format. One example can be viewed in the test_run folder provided here. It follows the structure presented here and matches the ISOGG tree nomenclature (https://isogg.org/tree/) with the first column containing haplogroup names of the corresponding tree branch:

<img src="/test_run/images/Custom_tree_hg.png" alt="Input file style" width="700"/>

SNPs that cannot be separated ("equal") are divided by commas in the same branch (green). SNPs of downstream branches are presented in the row below with one additional indentation using a tab (orange). And SNPs from parallel branches are on separate, mutually exclusive branches (blue).

#### 3.1.3) Output files
#### 3.1.3.1) phyloimputed.csv
The outcome are the observed (**D**,**A**,**X**) and imputed (**d**,**a**) allelic states for the initially reported SNPs complemented with the SNPs in the phylogenetic tree and a rooting SNP (ROOT). 

<img src="/test_run/images/Output_partly.png" alt="Input file style" width="350"/>

#### 3.1.3.2) haplogroups.csv
PhyloImpute cross-references the allelic states of all observed SNPs with the SNP relationships in the phylogenetic tree to verify the accuracy of the phylogenetic tree and the sequencing data. 

Following the tree branch with the most SNPs that are present in the analyzed sequence ("maximum resolution tree branch"), the haplogroup of the sample will be predicted. The output is given in the file **haplogroups.csv**, where each row present a sample, the predicted haplogroup and three statistical values to assess the reliability of the prediction. The first value is the **confidence value**, which counts the different derived alleles of SNPs from the **maximum resolution tree** branch divided by the total variants in that branch. A low **confidence value** can be simply due to low sequence coverage and does not mean that the predicted haplogroup is wrong. 

The second value, the **penalty value 1**, presents the number of markers from the **maximum resolution tree** branch in the input data with ancestral alleles divided by the total variants in that branch. The observed ancestral alleles could be the consequence of backmutations, but a high penalty value could hint towards an incorrect haplogroup prediction. 

The third value is the **penalty value 2** which counts the number of markers in parallel branches from the **maximum resolution tree** branch that were observed in the derived allelic state divided by the total number of variants in the tree. Markers in parallel branches in the derived allelic state could be the consequence of recurrent mutations that are identical by state, rather than by descent. However, a high penalty value 2 could indicate that the sequencing data comprises a mixture of DNA from different individuals.

#### 3.1.3.3) conflicting_SNPs.csv
Markers causing either of the two penalty values are stored in the additional file **conflicting_SNPs.csv**: Markers causing **penalty value 1** and **penalty value 2** are stored with the comments "(ancestral allele inside main branch)" and "(derived allele inside parallel branch)", respectively. (In the main output file (phyloimputed.csv), PhyloImpute keeps the observed allelic states.)


#
#### 3.2) VCF file
```
python PhyloImpute.py -input_format vcf -input ./test_run/input_vcf/ -output ./output -tree Y_minimal -vcf_ref GRCh37 -vcf_chr NC_000024.9
```
**Parameters:**

**-input_format** vcf

**-input** path to the folder with vcf files 

**-output** path to the folder for the output files

**-tree** path to the available phylogenetic tree {Y_minimal} [is mutually exclusive with -customtree] 

**-customtree** path to custom phylogenetic tree [is mutually exclusive with -tree] 

**-vcf_ref** Reference genome used for alignment {GRCh37, GRCh38}

**-vcf_chr** Chromosome nomenclature for chromosome of interest in the vcf file, e.g. NC_000024.9 for GRCh37

**-vcf_dic** path to dictionary file for the markers in the custom phylogenetic tree [exclusive to parameter -customtree]

#### 3.2.1) VCF Input file

Alternatively the path to a folder containing all vcf files can be provided. 

#### 3.2.2) Phylogenetic tree
#### 3.2.2.1) Pre-processed phylogenetic tree 
Currently, twp pre-processed phylogenetic trees are available for the human Y chromosome: The general Minimal Y tree (doi:10.1002/humu.22468) and haplogroup specific NAMQY tree (https://doi.org/10.1155/2024/3046495 ; Unpublished):
```
python PhyloImpute.py -input_format vcf -input ./test_run/input_vcf/ -output ./output -tree Y_minimal -vcf_ref GRCh37 -vcf_chr NC_000024.9
python PhyloImpute.py -input_format vcf -input ./test_run/input_vcf/ -output ./output -tree NAMQY -vcf_ref GRCh37 -vcf_chr NC_000024.9
```

#### 3.2.2.2) Custom phylogenetic tree
Alternatively, custom phylogenetic trees can be provided:
```
python PhyloImpute.py -input_format vcf -input ./test_run/input_vcf/ -output ./output -vcf_ref GRCh37 -vcf_chr NC_000024.9 -customtree ./test_run/minimal_y_tree_hgs_custom_example.csv -vcf_dic ./Y_minimal_dic.csv

```

The custom phylogenetic tree need to be made available by the user in the tab-separated _.csv_ format. One example can be viewed in the test_run folder provided here. It follows the structure presented here and matches the ISOGG tree nomenclature (https://isogg.org/tree/) with the first column containing haplogroup names of the corresponding tree branch:

<img src="/test_run/images/Custom_tree_hg.png" alt="Input file style" width="700"/>

SNPs that cannot be separated ("equal") are divided by commas in the same branch (green). SNPs of downstream branches are presented in the row below with one additional indentation using a tab (orange). And SNPs from parallel branches are on separate, mutually exclusive branches (blue).

##### 3.2.2.1) custom tree marker dictionary 
When using a custom tree, the user needs to provide information on the marker name, position, alleles in a dictionary file accessed with the parameter **-vcf_dic**

It needs to follow the presented structure:

<img src="/test_run/images/dictionary_file.png" alt="Input file style" width="350"/>

With the column names "marker" for the marker names identical to the ones used in the custom tree, the "GRCh37" and/or "GRCh38" columns (depending on your vcf file), "Anc" and "Der" columns for the ancestral and derived alleles of the respective variants. And the "Hg" column with the haplogroups defined by the respective marker.

#### 3.2.3) Output files
Same as above, see 3.1.3) Output files


# PhyloImpute
Impute missing data on non-recombining DNA using SNP phylogeny  
<br>

### 1) About PhyloImpute
PhyloImpute complements SNP data on non-recombining DNA, such as the human Y chromosome, by leveraging the SNPs' phylogeny as reported in a phylogenetic tree.  
<br>

### 2) Installation
Operating system: Linux

Type in the shell:
```bash
git clone https://github.com/ZehraKoksal/PhyloImpute.git
cd PhyloImpute/
python PhyloImpute.py -h
```  
<br>

### 3) Algorithm and Commands
### 3.1) Imputation
PhyloImpute imputes missing data by assuming that the SNPs in a clade of the phylogenetic tree leading up to a SNP with a derived allele are derived as well. SNPs on parallel branches are expected to be ancestral.

PhyloImpute can be run with a pre-processed csv input file or vcf files:  
<br>

#### 3.1.1) CSV file
```bash
python PhyloImpute.py -input_format csv -input ./test_run/testdata.csv -output ./output -tree Y_minimal
```

**Parameters:**

- **-input_format** csv  
- **-input** path to the input file  
- **-tree** path to the available phylogenetic tree {Y_minimal, NAMQY, ISOGG_2020} [mutually exclusive with -customtree]  
- **-customtree** path to custom phylogenetic tree [mutually exclusive with -tree]  
- **-output** path to the existing folder for the output files  
<br>

#### 3.1.1.1) CSV Input file
The user is required to provide the path to the input file in tab-separated `.csv` format.

<img src="/test_run/images/input_csv.png" alt="Input file style" width="350"/>

Rows represent variants, columns represent individuals.  
The header row contains individual labels (blue), and the second column contains variant names (orange).  
The table contains observed allelic states (green): ancestral **A**, derived **D**, or missing **X**.  
<br>

#### 3.1.1.2) Phylogenetic tree

##### 3.1.1.2.1) Pre-processed phylogenetic tree
Currently, a pre-processed phylogenetic tree is available for the human Y chromosome:
```bash
python PhyloImpute.py -input_format csv -input ./test_run/testdata.csv -output ./output -tree Y_minimal
```

##### 3.1.1.2.2) Custom phylogenetic tree
Alternatively, custom phylogenetic trees can be provided:
```bash
python PhyloImpute.py -input_format csv -input ./test_run/testdata.csv -output ./output -customtree ./test_run/minimal_y_tree_hgs_custom_example.csv
```

The custom phylogenetic tree must be provided in tab-separated `.csv` format.  
It follows the ISOGG tree nomenclature and this structure:

<img src="/test_run/images/Custom_tree_hg.png" alt="Input file style" width="700"/>

- SNPs that cannot be separated ("equal") are comma-separated in the same branch (green).  
- SNPs of downstream branches are indented using a tab (orange).  
- SNPs from parallel branches are in separate rows (blue).  
<br>

#### 3.1.1.3) Output files

##### 3.1.1.3.1) phyloimputed.csv
This file contains observed (**D**, **A**, **X**) and imputed (**d**, **a**) allelic states for reported SNPs plus additional SNPs from the phylogenetic tree, including a ROOT SNP.

<img src="/test_run/images/Output_partly.png" alt="Output preview" width="350"/>  
<br>

##### 3.1.1.3.2) haplogroups.csv
PhyloImpute compares allelic states of all observed SNPs with the SNP relationships in the phylogenetic tree to verify the accuracy of the phylogenetic tree and the sequencing data.

It outputs:
- Predicted haplogroup: based on tree branch with the most SNPs that are present in the analyzed sequence ("main tree branch")
- **Confidence value**: proportion of derived alleles in main tree branch. Low value can be due to low sequence coverage.  
- **Penalty value 1**: proportion of ancestral alleles in main branch. The observed ancestral alleles could be the consequence of backmutations, but a high penalty value could hint towards an incorrect haplogroup prediction.  
- **Penalty value 2**: proportion of derived alleles in parallel branches. Markers in parallel branches in the derived allelic state could be the consequence of recurrent mutations that are identical by state, rather than by descent. However, a high penalty value 2 could indicate that the sequencing data comprises a mixture of DNA from different individuals.    
<br>

##### 3.1.1.3.3) conflicting_SNPs.csv
Contains SNPs that cause penalty values:
- "(ancestral allele inside main branch)" → penalty value 1  
- "(derived allele inside parallel branch)" → penalty value 2

In the main output file (phyloimputed.csv), PhyloImpute keeps the observed allelic states.
<br>

#### 3.1.2) VCF file
```bash
python PhyloImpute.py -input_format vcf -input ./test_run/input_vcf/ -output ./output -tree Y_minimal -vcf_ref GRCh37 -vcf_chr NC_000024.9
```

**Parameters:**

- **-input_format** vcf  
- **-input** path to the folder with VCF files  
- **-output** output folder path  
- **-tree** phylogenetic tree {Y_minimal, NAMQY, ISOGG_2020} [exclusive with -customtree]  
- **-customtree** path to custom phylogenetic tree [exclusive with -tree]  
- **-vcf_ref** reference genome used {GRCh37, GRCh38, T2T}  
- **-vcf_chr** chromosome ID in VCF (e.g., NC_000024.9 for GRCh37)  
- **-vcf_dic** dictionary file for custom tree markers (used with -customtree)  
<br>

#### 3.1.2.1) VCF Input file
Provide a folder containing all `.vcf` files.  
<br>

#### 3.1.2.2) Phylogenetic tree

##### 3.1.2.2.1) Pre-processed phylogenetic tree
Supported trees:
- Minimal Y tree (doi:10.1002/humu.22468)  
- NAMQY tree (https://doi.org/10.1155/2024/3046495)  
- ISOGG 2020 tree (https://isogg.org/tree/)

```bash
python PhyloImpute.py -input_format vcf -input ./test_run/input_vcf/ -output ./output -tree Y_minimal -vcf_ref GRCh37 -vcf_chr NC_000024.9
python PhyloImpute.py -input_format vcf -input ./test_run/input_vcf/ -output ./output -tree NAMQY -vcf_ref GRCh37 -vcf_chr NC_000024.9
python PhyloImpute.py -input_format vcf -input ./test_run/input_vcf/ -output ./output -tree ISOGG_2020 -vcf_ref GRCh37 -vcf_chr NC_000024.9
```

##### 3.1.2.2.2) Custom phylogenetic tree
```bash
python PhyloImpute.py -input_format vcf -input ./test_run/input_vcf/ -output ./output -vcf_ref GRCh37 -vcf_chr NC_000024.9 -customtree ./test_run/minimal_y_tree_hgs_custom_example.csv -vcf_dic ./Y_minimal_dic.csv
```

Structure of the tree in tab-separated .csv file:
<img src="/test_run/images/Custom_tree_hg.png" alt="Tree file format" width="700"/>  

SNPs that cannot be separated ("equal") are divided by commas in the same branch (green). SNPs of downstream branches are presented in the row below with one additional indentation using a tab (orange). And SNPs from parallel branches are on separate, mutually exclusive branches (blue). 
One example can be viewed in the test_run folder provided here. 
<br>

##### 3.1.2.2.3) Custom tree marker dictionary
Required when using `-customtree`. Structure:

<img src="/test_run/images/dictionary_file.png" alt="Dictionary file format" width="350"/>

Columns:
- **marker**: names used in tree  
- **GRCh37**, **GRCh38**, **T2T**: positions  
- **Anc**, **Der**: alleles  
- **Hg**: haplogroups  
<br>

#### 3.1.2.3) Output files
Same as section 3.1.1.3  
<br>

### 3.2) Allele frequency plots
The output file of PhyloImpute (_phyloimputed.csv) can further be used to generate an allele frequency map of derived alleles of selected SNPs.
<img src="/test_run/images/freqmap_allelefreq_map_SNPX_basic.png" alt="general_allele_frequency_map" width="1000"/>

For this, the following code can be run:
```bash
python PhyloImpute.py -freqmap -input ./sample_data.csv -output ./freqmap -f_snp SNPX -f_coordinates ./sample_coordinates_example.csv -continent 'South America' 'North America' -af_map png
```

| **Flag**             | **Description**                                                                                                                                                                                                                                     | **Required/Optional** |
|----------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------|
| `-freqmap`           | Define this to generate allele frequency maps.                                                                                                                                                                                                     | Required               |
| `-input`             | Specify path to PhyloImpute output file (`_phyloimputed.csv`) or any file in the same format (see image in 3.1.1.3.1).                                                                                                                             | Required               |
| `-output`            | Define output file name.                                                                                                                                                                                                                            | Required               |
| `-f_snp`             | Define name of SNP for allele frequency map.                                                                                                                                                                                                       | Required               |
| `-f_coordinates`     | Provide path to tab-separated CSV file defining coordinates of samples from the `-input` file. Format: 1st column = sample names, 2nd = latitude, 3rd = longitude. Example: `/test_run/sample_coordinates_example.csv`.                             | Required               |
| `-continent`         | Specify one or several continents to plot: `Oceania`, `Africa`, `North America`, `Asia`, `South America`, `Europe`. Use single quotes, e.g., `'South America'`.                                                                                   | Optional               |
| `-country`           | Specify one or more countries to plot. Names must match **Countries_list.csv**. Use single quotes, e.g., `'Ecuador'`.                                                                                                                              | Optional               |
| `-whole_world`       | Plot the entire world map instead of specific regions.                                                                                                                                                                                             | Optional               |
| `-af_map`            | Select output file format for the allele frequency map: `svg`, `pdf`, or `png`. Default is `svg`.                                                                                                                                                  | Optional               |



### 3.2.1 Tune interpolator:
PhyloImpute maximizes the available information on the allelic states of SNPs by first imputing missing alleles (see above) and then by interpolating the remaining information between sample points using a radial basis function (RBF). This interpolation can be tuned by changing a parameter (epsilon). The default value of epsilon is 2.3, and the higher this value the stronger the "smoothing" of the data.
<br><br>
<img src="/test_run/images/PI_smoothing.png" alt="general_allele_frequency_map" width="1000"/>

I recommend illustrating datapoints with ancestral (black dots) and derived alleles (white dots) with the parameters **-ancestral_coordinates** and **-derived_coordinates**. Sometimes close-by datapoints can have extreme differences in allele frequencies (e.g., due to sampling strategy), which can cause artifacts in the RBF. The artifacts are visible as high derived allele frequencies of the SNP in a region, where no derived alleles are (i.e., no white dots). In these cases, the user should reduce the smoothing factor (epsilon) incrementally (by defining **-smoothing [number]**) until the artifact disappears.

```bash
python PhyloImpute.py -freqmap -input ./sample_data.csv -output ./freqmap -f_snp SNPX -f_coordinates ./sample_coordinates.csv -color pink -derived_coordinates -ancestral_coordinates -continent 'South America' -af_map png -smoothing 2
```

### 3.2.2 Additionally, some parameters can be changed to customize the maps:
<img src="/test_run/images/PI_plot_contours.png" alt="general_allele_frequency_map" width="1000"/>

- **-color**: The color palette can be changed by specifying one of these colors: blue,orange,pink,red,green,yellow,purple,violet,grey [Default:blue]
- **-contour**: Adjust the number of different shades for the allele frequencies [Default:15]

```bash
python PhyloImpute.py -freqmap -input ./sample_data.csv -output ./freqmap -f_snp SNPX -f_coordinates ./sample_coordinates.csv -color pink -continent 'South America' 'North America' -af_map png
```
### 3.2.3 More examples:
```bash
python PhyloImpute.py -freqmap -input ./sample_data.csv -output ./freqmap -f_snp SNP1 -f_coordinates ./sample_coordinates.csv -color orange -country 'Ecuador' -af_map png
```
<img src="/test_run/images/snp1_ecu.png" alt="country_specific_allele_frequency_map" width="500"/>
<br><br>

#don't use commas in snp names
import argparse
import time
import shutil
import glob
import os, sys
import gzip
import pandas as pd
import math
import numpy as np
import subprocess
import warnings
from datetime import datetime


print("\n")
start_time = time.time()

# Generate a timestamped filename
timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
log_filename = f"phyloimpute_log_{timestamp}.txt"

#User input paths
parser = argparse.ArgumentParser()

parser.add_argument("-input_format", choices=["vcf","csv"],help="Option 'vcf' is used when the available input format are .vcf or vcf.gz files. The vcf files need to be quality filtered already! The 'csv' option is used for tab-delimited csv file of a dataframe with individuals as columns and variants as rows. "
  "Insert variant data as 'A' for ancestral, 'D' for derived allele, and 'X' for missing data. Add variant names as first column and sequence name as header row. "
  "For both options: Refrain from using comma separated marker or sequence names.")

parser.add_argument("-input", required=True, help="For 'vcf' give the folder to the vcf files. For 'csv', give the path to one csv file. If -freqmap is specified, this needs to specify phyloimputed (or external) tab-separated csv file for allele frequency map.")
parser.add_argument('-v', action='store_true', help='Specify to output phyloimputed multi-sample vcf file.')
parser.add_argument("-output", required=True, help="Path to output folder.")
parser.add_argument("-vcf_ref", choices=["GRCh37","GRCh38","T2T"],help="When using the 'vcf' option, give the version of the reference genome used for the alignment: GRCh37, GRCh38, or T2T") 
parser.add_argument("-vcf_chr", help="How is the chromosome you want to analyze, named in your vcf file? This depends on the reference genome used and can be e.g. 'Y','chrY', 'NC_000024.9' for the human Y chromosome.")
parser.add_argument("-vcf_dic", help="If you are using a custom tree, you need to provide the path to a tab-separated .csv dictionary file that gives information on the genetic markers of interest that are also provided in the custom tree. Using the following column names: 'marker'(=marker names that are identical with the custom tree marker names),'GRCh37' and/or 'GRCh38' and/or 'T2T','Anc'(=ancestral allele), 'Der'(=derived allele)")
parser.add_argument("-tree", choices=["Y_minimal","NAMQY","ISOGG_2020"],help="Optional: path to tab-separated custom file of the phylogenetic SNP tree.")
parser.add_argument("-customtree", help="The path to a custom tree has to be provided in the same format as the tree files Y_minimal.csv, etc. Custom tree need to start with a root snp, e.g. 'ROOT'.")
parser.add_argument('-nucleotide', action='store_true', help='If specified, prints nucleotides (A,C,G,T) instead of allelic states (D/d,A/a,X) in phyloimputed file.')


parser.add_argument('-freqmap', action='store_true', help='Generate allele frequency map for specified SNP.')
parser.add_argument("-f_snp", help="Define SNP name to generate allele frequency maps.")
parser.add_argument("-f_coordinates", help="Provide path to csv file with sample names in first column (case-sensitive to input file) and coordinates in second column.")
parser.add_argument("-color", default="blue", choices=["blue","orange","pink","red","green","yellow","purple","violet","grey"],help="Select color [default:blue]") 
parser.add_argument("-contour", default="15", type=str, help="Define resolution of allele frequencies. Default: 15")
parser.add_argument("-smoothing", default="2.3", type=str, help="Define interpolated missing data points between sampling coordinates. The higher the value, the more extreme the smoothing. Default: 2.3")
parser.add_argument("-ancestral_coordinates", action="store_true",help="Specify if marking coordinates with ancestral alleles of SNP in black circles.")
parser.add_argument("-derived_coordinates", action="store_true",help="Specify if marking coordinates with ancestral alleles of SNP in white circles.")
parser.add_argument("-whole_world", action="store_true",help="Specify if visualizing map of whole world.")
parser.add_argument("-continent",choices=["Oceania","Africa","North America","Asia","South America","Europe"],nargs='+',help="Select one or more continents that you want to present. Add ' around each continent entry, e.g. 'South America'")
parser.add_argument("-country",nargs='+',help="Select one or more countries that you want to present. Nomenclature must align with the countries' nomenclature in file 'Countries_list.csv'.")
parser.add_argument("-af_map", default="svg", choices=["pdf","png","svg"],help="Specify output format of allele frequency map. Default: svg")


args = parser.parse_args()

current_dir = os.path.dirname(os.path.abspath(__file__))

#FOR VCF FORMAT OUTPUT FILES
def merge_vcfs(vcf_folder, output_vcf="merged.vcf.gz"):
    vcf_files = sorted(glob.glob(os.path.join(vcf_folder, "*.vcf"))) 
    if not vcf_files:
        raise FileNotFoundError("No VCF files found in the folder.")
    compressed_files = []
    for vcf in vcf_files:
        if vcf.endswith(".vcf"):
            gz = vcf + ".gz"
            with open(gz, "wb") as out:# Compress with bgzip, writing to gz
                subprocess.run(["bgzip", "-c", vcf], stdout=out, check=True)
            subprocess.run(["tabix", "-p", "vcf", gz], check=True)# Index with tabix
            compressed_files.append(gz)
    cmd = ["bcftools-1.20-rocky", "merge", "-Oz", "--force-samples", "-o", output_vcf] + compressed_files #"--force-samples", because of duplicate sample names
    subprocess.run(cmd, check=True)
    # Clean up
    to_remove = [f + ".gz" for f in vcf_files] + [f + ".gz.tbi" for f in vcf_files]
    subprocess.run(["rm", "-f"] + to_remove, check=True)
    print(f"Merged {len(compressed_files)} files into {output_vcf}")
    

#GENERATE ALLELE FREQUENCY MAP
if args.freqmap:
    if args.input_format:
        print("Warning: No need to provide 'input format'. Programme uses phyloimputed file as input data. This is not an error message, just information for user.")
    if not type(args.f_coordinates) == str:
        print("Please define sample names with corresponding coordinates using -f_coordinates. Provide path to csv file with sample names in first column (case-sensitive to input file) and coordinates in second column.")
        sys.exit()
    if not type(args.f_snp) == str:
        print("Please define SNP name to generate allele frequency map using -f_snp. Make sure to use the exact naming (case-sensitive) as in the input file.")
        sys.exit()
    cmd= ['python', 'freqmap.py',
            '-input', args.input,
            '-output', args.output,
            '-f_snp', args.f_snp,
            '-f_coordinates', args.f_coordinates,
            '-color', args.color,
            '-contour', args.contour,
            '-smoothing', args.smoothing,
            '-af_map', args.af_map
    ]
    if args.ancestral_coordinates:
        cmd.append('-ancestral_coordinates')    
    if args.derived_coordinates:
        cmd.append('-derived_coordinates')  
    if args.whole_world:
        cmd.append('-whole_world')  
    elif args.continent:
        cmd.extend(['-continent'] + args.continent) 
    elif args.country:
        cmd.extend(['-country'] + args.country) 
    subprocess.run(cmd)    
    sys.exit()


if not args.input_format:
    print("Please define '-input_format' to run PhyloImpute.")
    sys.exit()
    
#FOR VCF INPUT FILE
print("PhyloImpute v1.3 is running...")
print("\n")
with open(log_filename, "w") as f:
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    f.write(f"[{timestamp}]\tPhyloImpute v1.3 started ...\n")
    
def check_for_vcf_files(directory):
    vcf_files = glob.glob(os.path.join(directory, '*.vcf'))
    if vcf_files:
        print(f'Found vcf files.')
        with open(log_filename, "a") as f:
            timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
            f.write(f"[{timestamp}]\tFound .vcf files: {vcf_files}\n")
        return vcf_files
    else:
        print('\nNo .vcf files found.')
        return []

#Set this as default to avoid errors
if args.input_format == "csv":
    args.vcf_ref="GRCh37"

if args.tree=="Y_minimal":
    dic= pd.read_csv("./Y_minimal_dic.csv", sep="\t")
    #list of positions frm dictionary to filter from vcf file
    dic[args.vcf_ref] = dic[args.vcf_ref].astype('Int64')
    dic[args.vcf_ref] = dic[args.vcf_ref].astype(str)
    dic_pos = dic[args.vcf_ref].values.tolist()
    dic_Anc = dic["Anc"].values.tolist()
    dic_Der = dic["Der"].values.tolist()
    dic_marker = dic["marker"].values.tolist()
elif args.tree=="NAMQY":
    dic= pd.read_csv("./NAMQY_dic.csv", sep="\t")
    #list of positions frm dictionary to filter from vcf file
    dic[args.vcf_ref] = dic[args.vcf_ref].astype('Int64')
    dic[args.vcf_ref] = dic[args.vcf_ref].astype(str)
    dic_pos = dic[args.vcf_ref].values.tolist()
    dic_Anc = dic["Anc"].values.tolist()
    dic_Der = dic["Der"].values.tolist()
    dic_marker = dic["marker"].values.tolist()
elif args.tree=="ISOGG_2020":
    dic= pd.read_csv("./ISOGG_2020_dic.csv", sep="\t")
    #list of positions frm dictionary to filter from vcf file
    dic[args.vcf_ref] = dic[args.vcf_ref].astype('Int64')
    dic[args.vcf_ref] = dic[args.vcf_ref].astype(str)
    dic_pos = dic[args.vcf_ref].values.tolist()
    dic_Anc = dic["Anc"].values.tolist()
    dic_Der = dic["Der"].values.tolist()
    dic_marker = dic["marker"].values.tolist()
#If you are using a custom tree you need to also provide a custom dictionary and extract the information here
elif type(args.customtree)==str:
    dic= pd.read_csv(args.vcf_dic, sep="\t")
    #list of positions frm dictionary to filter from vcf file
    dic[args.vcf_ref] = dic[args.vcf_ref].astype('Int64')
    dic[args.vcf_ref] = dic[args.vcf_ref].astype(str)
    dic_pos = dic[args.vcf_ref].values.tolist()
    dic_Anc = dic["Anc"].values.tolist()
    dic_Der = dic["Der"].values.tolist()
    dic_marker = dic["marker"].values.tolist()


#Make dic for later with SNP name and hg
marker_to_hg = dict(zip(dic["marker"], dic["Hg"]))    
     
        
#get_marker_info_from_vcf by using a dictionary file
if args.input_format == "vcf":
    
        
    #Separate dic file from ref
    ref_dic = dic.copy()
    dic = dic.drop(columns=["Ref"])
    

    #Make input table for Phyloimpute
    #turn positions to strings, so they can be recognized later in the data tables
    input_table = pd.DataFrame({'pos':dic_pos, 'Sample':dic_marker, 'der': dic_Der})
    input_table_raw = pd.DataFrame({'pos':dic_pos, 'Sample':dic_marker, 'der': dic_Der})
    input_table = input_table.drop_duplicates()
    input_table = input_table.sort_index().reset_index(drop=True)
    column_missing_data = len(input_table)*["X"]
    # print(input_table)
    # print(ee)


def extract_sample_name(file_path):
    split_string = file_path.split("/")
    sample_name_format = split_string[-1]
    sample_name = sample_name_format.replace(".csv", "")
    sample_name = sample_name_format.replace(".vcf", "")
    return sample_name


def filter_vcf_file(file_path):
    data = []
    columns = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#CHROM'):
                #Split header line to get column names
                columns = line.strip().split('\t')
            elif not line.startswith('##'):
                data.append(line.strip().split('\t'))
    #Create dataframe from the data
    df = pd.DataFrame(data,columns=columns)

    # df = df.drop_duplicates(subset=["POS"])
    # print(df)
    #Filter chromosome column
    df = df[df["#CHROM"].isin([args.vcf_chr])]
    
    #Keep only first base in ref and of alts
    df['REF'] = df['REF'].str[0]
    # Split by comma
    df['ALT'] = df['ALT'].str.split(',')

    # Keep only first letter per allele
    df['ALT'] = df['ALT'].apply(lambda alleles: [a[0] for a in alleles])

    #Filter positions of variants from database
    df = df[df["POS"].isin(dic_pos)]
    

    #check if multi-sample vcf
    # if len(df.columns)>10:
        # print("Multi-sample vcf recognized.")
    samples_list = df.columns[9:].tolist()
    print(samples_list)
    # print(Ee)
    

    #Also filter that is the derived allele that we see:
    ##Make list of variant positions in the same order
    vcf_pos_list = df["POS"].tolist()
    # print(vcf_pos_list)
    dic_vcf_pos = dic[dic[args.vcf_ref].isin(vcf_pos_list)]
    # print(dic_vcf_pos)
    # print(ee)

    warnings.simplefilter(action='ignore', category=pd.errors.SettingWithCopyWarning)
    #now sort according to order in the list
    dic_vcf_pos.loc[:, 'sort_key']= pd.Categorical(dic_vcf_pos.loc[:, args.vcf_ref], categories=vcf_pos_list, ordered=True)
    
    
    sorted_dic_vcf_pos= dic_vcf_pos.sort_values('sort_key').drop('sort_key', axis=1)

    
    #and remove index
    sorted_dic_vcf_pos.reset_index(drop=True, inplace=True)
    #Take the derived alleles from dictionary in same order to the vcf file 
    #and remove index of filtered vcf
    df.reset_index(drop=True, inplace=True)
    #add derived alleles from dic per marker
    sorted_dic_vcf_pos_2 = pd.DataFrame()
    # sorted_dic_vcf_pos_2[["POS","der_from_dic_for_pos"]] = sorted_dic_vcf_pos[[args.vcf_ref,"Der"]]
    sorted_dic_vcf_pos_2[["POS","der_from_dic_for_pos","anc_dic"]] = sorted_dic_vcf_pos[[args.vcf_ref,"Der","Anc"]]
    sorted_dic_vcf_pos_2 = sorted_dic_vcf_pos_2.drop_duplicates()

    df = pd.merge(df,sorted_dic_vcf_pos_2, on="POS", how="left")

    # print(df)
    # print(df.iloc[:, 9])
    # print(ee)
    
    #for multi-sample vcfs make dfs to store information in loop
    df_filtered_list_D_all = []
    df_filtered_list_D_der_all = []
    df_filtered_list_A_all = []
    df_filtered_list_A_der_all = []

    for sample, snr in zip(samples_list, range(len(samples_list))):
        print(sample)
        # print(snr)
        column = 9 + snr
        # print(column)
        # for vcfs that also contain non polymorphic sites: if alt allele matches derived allele AND genotype column (this is index 9 + additional samples (for multisample vcf) 
        df_filtered_1 = df[
        (df['ALT'].str[0] == df['der_from_dic_for_pos']) &
        (df.iloc[:, column].astype(str).str.startswith("1"))]
        
        df_filtered_0 = df[
        (df['ALT'].str[1] == df['der_from_dic_for_pos']) &
        (df.iloc[:, column].astype(str).str.startswith("2"))]
        
        df_filtered_2 = df[
        (df['REF'] == df['der_from_dic_for_pos']) &
        (df.iloc[:, column].astype(str).str.startswith("0"))]
        df_filtered = pd.concat([df_filtered_1, df_filtered_2, df_filtered_0], ignore_index=True)  

        # print(df_filtered)

        #Make list of the position column
        df_filtered_list_D = df_filtered["POS"].values.tolist()
        df_filtered_list_D_der = df_filtered["der_from_dic_for_pos"].values.tolist()

        
        df_filtered_A1 = df[
        (df['ALT'].str[0] == df['anc_dic']) &
        (df.iloc[:, column].astype(str).str.startswith("1"))]
        
        df_filtered_A0 = df[
        (df['ALT'].str[1] == df['anc_dic']) &
        (df.iloc[:, column].astype(str).str.startswith("2"))]
        
        df_filtered_A2 = df[
        (df['REF'] == df['anc_dic']) &
        (df.iloc[:, column].astype(str).str.startswith("0"))]
        df_filtered_A = pd.concat([df_filtered_A1, df_filtered_A2, df_filtered_A0], ignore_index=True)
            
  
        #Make list of the position column
        df_filtered_list_A = df_filtered_A["POS"].values.tolist()
        df_filtered_list_A_der = df_filtered_A["der_from_dic_for_pos"].values.tolist()
        
        df_filtered_list_known = df_filtered_list_A + df_filtered_list_D
        
        # print(len(df_filtered_list_D))
        # print(len(df_filtered_list_A))
        # print(len(df_filtered_list_known))
        
        #for now: ignore
        # df_unknown = df[~df["POS"].isin(df_filtered_list_known)]
        # print(df_unknown.iloc[:,3:])
        # df_unknown = df_unknown[~df_unknown[sample].str.startswith(".")]
        # print(df_unknown.iloc[:,3:10])


        df_filtered_list_D_all.append(df_filtered_list_D)
        df_filtered_list_D_der_all.append(df_filtered_list_D_der)
        df_filtered_list_A_all.append(df_filtered_list_A)
        df_filtered_list_A_der_all.append(df_filtered_list_A_der)
        
        # print(df_filtered_list_D_all)
        # print(eee)

    return df_filtered_list_D_all, df_filtered_list_D_der_all, df_filtered_list_A_all, df_filtered_list_A_der_all, samples_list



if args.input_format == "vcf":
    vcf_files = check_for_vcf_files(args.input)
    if type(args.vcf_ref)==str and type(args.vcf_chr)==str:
        print("\nNecessary input files are provided. Processing files now!\n")
        #Now that checking is over, we loop through 
        for vcf_file in vcf_files:
            with open(log_filename, "a") as f:
                timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
                f.write(f'[{timestamp}]\tProcessing file: {vcf_file}\n')
            derived_variant_pos_all, derived_variant_der_all, ancestral_variant_pos_all, ancestral_variant_der_all, samples_list =filter_vcf_file(vcf_file)
            # print(derived_variant_pos_all)
            
            # print(samples_list)
            input_table['pos_der'] = list(zip(input_table['pos'], input_table['der']))
            # print(input_table)
            # print(ee)
            for sample_name, derived_variant_pos, derived_variant_der, ancestral_variant_pos, ancestral_variant_der in zip(samples_list, derived_variant_pos_all, derived_variant_der_all, ancestral_variant_pos_all, ancestral_variant_der_all):
                # sample_name = extract_sample_name(vcf_file)
                #Add sample column to input table with "X"s in it
                input_table[sample_name] = column_missing_data
                # print(input_table)
                remove_set_d = set(zip(derived_variant_pos, derived_variant_der))
                remove_set_a = set(zip(ancestral_variant_pos, ancestral_variant_der))
                
                # Assign values directly using .loc and isin
                mask_d = input_table['pos_der'].isin(remove_set_d)
                mask_a = input_table['pos_der'].isin(remove_set_a)

                input_table.loc[mask_d, sample_name] = "D"
                input_table.loc[mask_a, sample_name] = "A"
                
                # input_table_D = input_table[input_table['pos_der'].isin(remove_set_d)]
                # input_table_D.loc[:,sample_name]="D"
                # input_table_A = input_table[input_table['pos_der'].isin(remove_set_a)]
                # input_table_A.loc[:,sample_name]="A"
                #Combine dfs
                # input_table_no_info = input_table[~input_table['pos_der'].isin(remove_set_d)]
                #Repeat for ancestral
                # input_table_no_info = input_table_no_info[~input_table_no_info['pos_der'].isin(remove_set_a)]
                # input_table = pd.concat([input_table_D, input_table_A, input_table_no_info])
                # print(input_table)

            
            #clean up            
            input_table = input_table.drop(columns='pos_der')
            input_table = input_table.drop_duplicates(subset=["Sample","pos"])
            input_table = input_table.sort_index().reset_index(drop=True)
            
            # print(input_table)
            # print(ee)

        #replace pos column and reindex
        input_table.reset_index(drop=True, inplace=True)
        input_table["pos"]=range(1,len(input_table)+1)
        df = input_table
        # print(df)
        df.set_index("Sample", inplace=True)
        #rename column pos
        df.rename(columns={'pos':0}, inplace=True)
        df = df.drop(columns=['der'])
        # print(df)
        # print(ee)
    else:
        print("\nOne or both of the input files vcf_ref or vcf_chr are missing! Please provide and rerun. See help function for more information. \n")
        with open(log_filename, "a") as f:
            timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
            f.write(f'[{timestamp}]\tOne or both of the input files vcf_ref or vcf_chr are missing! Please provide and rerun. See help function for more information.\n')


if args.input_format == "csv":
    #import sample dataframe and variant order 
    df = pd.read_csv(args.input, sep="\t", index_col=1) #INPUT CHANGE
    


#WITH INPUT TABLE

#skip numbering column
df = df.iloc[:,1:]
#replace missing data with X and (*) with X
df = df.fillna("X")
df = df.replace(to_replace=r'\(.*\)', value='X',regex=True)

#TREE to DATABASE
if type(args.customtree) == str:
    print("CUSTOM TREE")
    print(args.customtree)
    tree_og = pd.read_csv(args.customtree, sep="\t", header=None)
    hgs = tree_og.iloc[1:,0]
    hgs_list = hgs.values.tolist()
    tree_df = tree_og.iloc[:,1:]
    tree_df.columns = range(len(tree_df.columns)) #reindex columns starting from 0
    for col_name, col in tree_df.items():
        first_non_nan_index = col.first_valid_index() ##Get the index with the first non nan
        for index, value in col.iloc[first_non_nan_index:].items():#tree.iterrows():
            if index==len(tree_df)-1:
                break
            if pd.notna(index):
                marker=value
            if col_name==0:
                tree_df.loc[:,col_name]  = col.fillna(marker) #replace all nan with the first marker non-nan
                break #so we can skip to the next column
            elif pd.isna(col.iloc[index+1]): 
                prev_upstream=tree_df.loc[index,col_name-1]
                if tree_df.loc[index+1,col_name-1]==prev_upstream:
                    warnings.simplefilter(action='ignore', category=pd.errors.SettingWithCopyWarning)
                    tree_df.loc[index+1,col_name]  = marker
                    if index==len(tree_df)-2:
                        break

    tree_df= tree_df.iloc[1:,:] #skip the row with only "ROOT" entry
    tree_df.insert(loc=0, column="hgs", value=hgs_list)
    custom_tree_path = current_dir + "/custom_tree.csv"
    tree_df.to_csv(custom_tree_path, index=False, sep='\t')
    print("SNP database generated and saved!")
    with open(log_filename, "a") as f:
        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        f.write(f'[{timestamp}]\tSNP database for custom tree generated and saved!\n')


if args.tree == "ISOGG_2020":
    filename = args.tree + ".csv.gz"
    tree_csv_path = os.path.join(current_dir,filename)
    tree_df = pd.read_csv(tree_csv_path, sep="\t", compression='infer')
elif type(args.tree) == str:
    filename = args.tree + ".csv"
    tree_csv_path = os.path.join(current_dir,filename)
    tree_df = pd.read_csv(tree_csv_path, sep="\t")
    #filter out the hg column 
tree= tree_df.iloc[:,1:]
    #now get the haplogroup column
hg_table = tree_df.iloc[:,0]
hg_table = hg_table.values.tolist()

#turn rows into lists
tree_lists=[row.dropna().tolist() for index, row in tree.iterrows()]


#get all snps in the whole tree
tree_snps = []
for i, r in tree.iterrows():
    for e in r:
        if ',' in str(e):
            tree_snps.extend(e.strip() for e in str(e).split(',') if e.strip())
        else:
            tree_snps.append(e)
tree_snps = set(tree_snps) 
tree_snps = {item for item in tree_snps if not (isinstance(item, float) and math.isnan(item))}

# Add all snp names that are not in the sample df, as rows filled with X's
df_snps = set(df.index.tolist())

adding_snps = tree_snps - df_snps
adding_snps = {x for x in adding_snps if x is not np.nan}
adding_snps = {x for x in adding_snps if x != ""}
df_exclusive_SNPs = df_snps - tree_snps


#add the missing snps from the tree to the data
df_adding_snps = pd.DataFrame(index=list(adding_snps), columns=df.columns, data="X")
df = pd.concat([df,df_adding_snps])

# Temporarily move index to a column to remove duplicate entries of the same marker
df = df.reset_index().drop_duplicates()
df = df.set_index("index")
df.index.name = None


#B) START REPLACING SAMPLE DATAFRAME
def preprocess_set(s):
    processed_set= set()
    for item in s:
        if "," in item:
            components=item.split(",")
            processed_set.update(components)
        else:
            processed_set.add(item)
    return processed_set
questionable_SNPs= [] #something is off with these variants position

Haplogroup_info = ["Sample\tHg\tConfidence_value\tPenalty_value1\tPenalty_value2\tdownstream_ancestral\tdownstream_unknown"]
for col_name, col in df.items(): #iteritems()
    lengths_list = []
    # print(df[col_name])
    filtered_D= df[col_name]=="D"    
    sample_D_list = df[filtered_D].index.tolist()#make a list of the index values (SNP names) of the filtered dataframe
    # print(filtered_D)
    # print(sample_D_list)
    # print(ee)
    
    if len(sample_D_list) == 0:
        print(f"No derived SNPs were found for {col_name}, therefore no imputation is possible. You could try using a different phylogenetic tree for imputation or retrieve more sequencing data.\n")
        with open(log_filename, "a") as f:
            timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
            f.write(f'[{timestamp}]\tNo derived SNPs were found for {col_name}, therefore no imputation is possible. You could try using a different phylogenetic tree for imputation or retrieve more sequencing data.")\n')
        #FOR HG
        Hg_sample = "No hg pred possible, because no informative derived alleles"
        Hg_support = "-"
        Hg_penalty1 = "-"
        Hg_penalty2 = "-"
        dwn_A = "-"
        dwn_X = "-"
        Haplogroup_info.append(col_name+"\t"+Hg_sample+"\t"+str(Hg_support)+"\t"+str(Hg_penalty1)+"\t"+str(Hg_penalty2)+"\t"+str(dwn_A)+"\t"+str(dwn_X)) #col_name is the sample name
        
        #For questionable SNPs
        questionable_entry_A = str(col_name)+": - "
        questionable_SNPs.append(questionable_entry_A)
    else:
        filtered_A= df[col_name]=="A"
        sample_A_list = df[filtered_A].index.tolist()#make a list of the index values (SNP names) of the filtered dataframe
        sample_A_set=set(sample_A_list)
        # print(sample_A_set)
        #Compare each sample's derived lists with the database lists
        # print(Ee)
        processed_D_list=preprocess_set(sample_D_list)
        for db_list in tree_lists:
            processed_db_list=preprocess_set(db_list)
            overlapping_entries=processed_D_list.intersection(processed_db_list)
            len_overlap = len(overlapping_entries)
            lengths_list.append(len_overlap)
        max_index = lengths_list.index(max(lengths_list))
        
        #get the database row that will be used to replace missing data
        #are any of these items overlapping with ancestral observations in the sample?
        new_D_db =set(tree_lists[max_index])  # matches = {item for item in new_D_db if "F313" in item}
        # print(new_D_db)
        highest_resolution_D_tree = tree_lists[max_index][-1]
        highest_resolution_D_tree = set(highest_resolution_D_tree.split(","))
        # print(highest_resolution_D_tree)
        
        highest_resolution_SNP=processed_D_list.intersection(highest_resolution_D_tree)

        processed_set1=preprocess_set(sample_A_list)
        max_d=preprocess_set(tree_lists[max_index])

        
        #GET INFO ON HG OF SAMPLE           ##HG
        overlapping_derived=processed_D_list.intersection(max_d) #overlap between derived snps and snps of max branches
        Hg_sample = hg_table[max_index]
        
        #Now check if downstream markers are ancestral or absent

        if len(overlapping_derived) == 0:
            print(f"No derived SNPs in the phylogenetic tree were found for {col_name}, therefore no imputation is possible. You could try using a different phylogenetic tree for imputation or retrieve more sequencing data.\n")
            with open(log_filename, "a") as f:
                timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
                f.write(f'[{timestamp}]\tNo derived SNPs were found for {col_name}, therefore no imputation is possible. You could try using a different phylogenetic tree for imputation or retrieve more sequencing data.")\n')
            #FOR HG
            Hg_sample = "No hg pred possible, because no informative derived alleles"
            Hg_support = "-"
            Hg_penalty1 = "-"
            Hg_penalty2 = "-"
            dwn_A = "-"
            dwn_X = "-"
            Haplogroup_info.append(col_name+"\t"+Hg_sample+"\t"+str(Hg_support)+"\t"+str(Hg_penalty1)+"\t"+str(Hg_penalty2)+"\t"+str(dwn_A)+"\t"+str(dwn_X)) #col_name is the sample name
            
            #For questionable SNPs
            questionable_entry_A = str(col_name)+": - "
            questionable_SNPs.append(questionable_entry_A)
        
        else:
        
            # matches = [lst for lst in tree_lists if all(x in lst for x in new_D_db)]
            matches = [lst for lst, lst_set in zip(tree_lists, map(set, tree_lists))
            if new_D_db.issubset(lst_set) and lst_set != new_D_db]

            
            #removed all entries up to most resolved SNP, giving downstream SNPs
            filtered_matches = [
                [x for x in lst if x not in new_D_db]
                for lst in matches
            ]

            
            flat_list_dwnstream = [x for lst in filtered_matches for x in lst]

            
            #remove duplicates
            flat_list_dwnstream = list(set(flat_list_dwnstream))
            # print(flat_list_dwnstream)
            #split by commas
            flat_list_dwnstream_preprocess = preprocess_set(flat_list_dwnstream)

            dwnstream_ancestral = processed_set1.intersection(flat_list_dwnstream_preprocess)
            dwnstream_ancestral = {x for x in dwnstream_ancestral if x != ""}

            dwnstream_X = flat_list_dwnstream_preprocess - dwnstream_ancestral
            dwnstream_X = {x for x in dwnstream_X if x != ""}

            
            #Look for haplogroups of these SNPs
            Hgs_Anc_down = [marker_to_hg[m] for m in dwnstream_ancestral if m in marker_to_hg]
            Hgs_Anc_down_trimmed = [item.replace(Hg_sample, "~") for item in Hgs_Anc_down]

            
            #Combine Hgs and SNPs
            dwn_A = [f"{a}-{b}" for a, b in zip(Hgs_Anc_down_trimmed, dwnstream_ancestral)]

            
            #Repeat for missing SNPs
            Hgs_X_down = [marker_to_hg[m] for m in dwnstream_X if m in marker_to_hg]
            Hgs_X_down_trimmed = [item.replace(Hg_sample, "~") for item in Hgs_X_down]
            dwn_X = [f"{a}-{b}" for a, b in zip(Hgs_X_down_trimmed, dwnstream_X)]
            # print(dwn_X) #OUTPUT
            

            
            Hg_support = str(len(overlapping_derived))+"/"+str(len(max_d))
            #Penalty comprises of variants in the max_d_branch that are ancestral (but are expected to be derived)
            #But also of variants that are derived but not in the max_d_branch
            Hg_penalty_part1 = processed_set1.intersection(max_d)#How many are ancestral in our max res branch
            Hg_penalty1 = str(len(Hg_penalty_part1))+"/"+str(len(max_d))#How many are ancestral
            
            Hg_penalty_part2 = processed_D_list - max_d #How many in other branches are derived #TO DO: exclude SNPs that are not in tree
            Hg_penalty_part2 = Hg_penalty_part2 - df_exclusive_SNPs #remove SNPs that are derived in data, but not in tree
            
            #uncertain
            if int(len(Hg_penalty_part1)) > int(len(overlapping_derived)):
                Hg_sample = "*" + Hg_sample

            Hg_sample = Hg_sample + "-" + "/".join(highest_resolution_SNP)
            

            if "ROOT" in tree_snps:
                total_snps_tree = len(tree_snps)-2 #Remove title and ROOT
            else:
                total_snps_tree = len(tree_snps)-1 #Remove title, No ROOT
            Hg_penalty2 = str(len(Hg_penalty_part2))+"/"+str(total_snps_tree) #divided by total number of variants
            
            # Haplogroup_info.append(col_name+"\t"+Hg_sample+"\t"+str(Hg_support)+"\t"+str(Hg_penalty1)+"\t"+str(Hg_penalty2)) #col_name is the sample name
            Haplogroup_info.append(col_name+"\t"+Hg_sample+"\t"+str(Hg_support)+"\t"+str(Hg_penalty1)+"\t"+str(Hg_penalty2)+"\t"+str(dwn_A)+"\t"+str(dwn_X))
        
            # overlapping_entries=processed_set1.intersection(max_d)
            if len(processed_D_list) > 0:
                if len(Hg_penalty_part1)>0:
                    questionable_entry_A = str(col_name)+":"+str(Hg_penalty_part1)+" (ancestral allele inside main branch)"
                    questionable_SNPs.append(questionable_entry_A)
                if len(Hg_penalty_part2)>0:
                    questionable_entry_D = str(col_name)+":"+str(Hg_penalty_part2)+" (derived allele inside parallel branch)"
                    questionable_SNPs.append(questionable_entry_D)

            
        #A)get the SNPs from the database row, but exclude those that were ancestral in the sample
        to_replace_derived = max_d.difference(processed_set1)
        new_D_db = max_d-processed_set1-processed_D_list #the snps in the database row excluding the ancestral and derived variants from the current sample
        
        if len(new_D_db) > 0:
            filtered_df = df[df.index.isin(new_D_db)]
            unfiltered_df = df[~df.index.isin(new_D_db)] #ALL OTHER VARIANTS to combine in the end with the altered filtered_df 
            filtered_df.loc[:,col_name]="d" #make the all derived (filling gaps step) #*for inference
            df=pd.concat([filtered_df,unfiltered_df])
        
        #B)turn X to A
        A_list = []
        for i in range(len(tree_lists[max_index])):
            #filter for each column between max_d list and dataframe tree
            df_mismatch = tree[~tree.iloc[:,i].isin([tree_lists[max_index][i]])]
            #if there are parallel branches, df_mismatch will have a length > 0
            if len(df_mismatch) > 0:
                #take all variants that are parallel to our sample based on the one column
                data_array=df_mismatch.values
                parallel_markers=data_array.flatten().tolist()
                #remove nan
                parallel_markers = {x for x in parallel_markers if pd.notnull(x)}
                parallel_markers = {x for x in parallel_markers if x is not np.nan}

                A_list.append(list(parallel_markers))
                #to set, to remove duplicate entries   
        A_list=set([item for sublist in A_list for item in sublist]) #flatten list in lists
        A_list = A_list - set(tree_lists[max_index]) #set(tree_lists[max_index]) are the derived snps
        A_list_markers=preprocess_set(A_list)
        #also remove the derived single (already processed snps) that are derived in the sample:
        A_list_markers = A_list_markers -set(sample_D_list)
        #now remove the ones that are already known to be ancestral
        A_list_markers=A_list_markers-sample_A_set
        #for recurrent snps: remove the ones that are imputed to be derived
        A_list_markers=A_list_markers-max_d

        if len(A_list_markers) > 0:
            filtered_df = df[df.index.isin(A_list_markers)]
            unfiltered_df = df[~df.index.isin(A_list_markers)] #ALL OTHER VARIANTS to combine in the end with the altered filtered_df 
            filtered_df.loc[:,col_name]="a" #make the all derived (filling gaps step) #*for inference
            df=pd.concat([filtered_df,unfiltered_df])

#flatten list of sets:
questionable_SNPs=list(questionable_SNPs)
#remove "snp" called ROOT (that is not a real snp)
df = df[df.index != "ROOT"]

def get_variant_code(row, samp):
    # print(row)
    val = row[f"{samp}_y"] #_y
    # print(val)
    # print(ee)
    if pd.isna(val):
        return val  # Keep NaN
    elif val == row["REF"]:
        return "0:999"
    else: #allele is alt
        alt_alleles = str(row["ALT"]).split(",")
        try:
            idx = alt_alleles.index(val)
            return f"{idx + 1}:999"  # ALT indices start at 1
        except ValueError:
            return val  # Not found in ALT list → leave unchanged

#If user wants nucleotides instead of allelic states
if args.v == True:
    vcf_input_folder = args.input
    output_path = args.output + "_phyloimputed.vcf.gz"
    with open(log_filename, "a") as f:
        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        f.write(f'[{timestamp}]\tStart merging vcf files\n')
    merge_vcfs(vcf_input_folder, output_path)
    with open(log_filename, "a") as f:
        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        f.write(f'[{timestamp}]\tSuccessfully merged vcf files\n')
    #Open as df
    data = []
    columns = []
    with gzip.open(output_path, 'rt') as file:  # 'rt' = read text mode
        for line in file:
            if line.startswith("##bcftools_mergeCommand=merge"):
                #To make sure the sample column headers match the vcf files (not the original BAM files)
                line = line.replace("##bcftools_mergeCommand=merge -Oz --force-samples -o", "")
                lst = line.split()
                lst = lst[1:-5]
                lst = [item.replace(".vcf.gz;","") for item in lst]
                lst = [item.replace(".vcf.gz","") for item in lst]
                lst = [item.split("/")[-1] for item in lst]
            elif line.startswith('#CHROM'):
                # Split header line to get column names
                columns = line.strip().split('\t')
                #replace sample names to match the vcf files
                # for i in range(len(lst)):
                    # i = i +1
                    # columns[-i] = lst[-i]
            elif not line.startswith('##'):
                data.append(line.strip().split('\t'))
    merged_vcf = pd.DataFrame(data, columns=columns)
    # print(output_path)
    # print(merged_vcf)
    # print(ee)
    
    merged_vcf_pos = merged_vcf.loc[merged_vcf['#CHROM'] == args.vcf_chr, 'POS'].tolist()
    # merged_vcf_pos = merged_vcf["POS"].values.tolist()
    # print(merged_vcf_pos)
    # print(ee)
    
    #remove unimputed alleles and missing alleles
    # pd.set_option('future.no_silent_downcasting', True)
    df_trimmed = df.replace({'A': np.nan, 'D': np.nan, 'X': np.nan})
    # print(df_trimmed)
    
    # print(ee)
    
    df_trimmed = df_trimmed.dropna(how='all') #remove rows that are completly NaN in all samples
    # print(df_trimmed)
    
    dic_relevant = pd.DataFrame()
    dic_relevant[[args.vcf_ref, "Anc", "Der","marker"]] = dic[[args.vcf_ref, "Anc", "Der","marker"]]
    dic_relevant = dic_relevant.drop_duplicates()
    # print(dic_relevant)
    # print(ee)
    
    
    #Merge trimmed df and dic_relevant
    df_dic_merge = dic_relevant.merge(df_trimmed, left_on="marker", right_index=True, how= "inner")
    # print(df_dic_merge)
    # print(Ee)
    
    for col in df_dic_merge.columns[4:]:  # from 5th column onward
        df_dic_merge[col] = df_dic_merge.apply(
            lambda row: row["Anc"] if row[col] == "a" else row["Der"] if row[col] == "d" else row[col],
            axis=1)
    # print(df_dic_merge)
    # print(Ee)
    df_dic_merge = df_dic_merge.drop(columns=["Anc","Der"])
   

    #split imputation file into positions that are already in vcf file and those that are not
    df_dic_merge_exi = df_dic_merge[df_dic_merge[args.vcf_ref].isin(merged_vcf_pos)]
    df_dic_merge_new = df_dic_merge[~df_dic_merge[args.vcf_ref].isin(merged_vcf_pos)]
    samp_list = list(df_dic_merge_exi.columns[2:])
    # print(df_dic_merge_exi)
    # print(df_dic_merge_new)
    # print(ee)
    

    #Process existing positions
    #merge vcf and existing imp snps
    merge_vcf_imp_exi = merged_vcf.merge(df_dic_merge_exi, left_on="POS", right_on=str(args.vcf_ref), how= "left")
    # print(merge_vcf_imp_exi)
    # print(Ee)
    
    #Now compare the info from original _x and imputed _y column of each sample/panel vcf
    
    #Update Alt column of merge_vcf_imp_exi
    samp_y_cols = [col for col in merge_vcf_imp_exi.columns if col.endswith('_y')]
    # print(samp_y_cols)
    # print(Ee)
    
    
    def fill_alt_from_y(row, samp_y_cols):
        if row["ALT"] != ".":
            return row["ALT"]  # Keep existing ALT

        # Gather all values from *_y columns that are not equal to REF and not NaN
        alternatives = {
            row[col]
            for col in samp_y_cols
            if pd.notna(row[col]) and row[col] != row["REF"]
        }

        if not alternatives:
            return "."  # Keep as-is if no valid alternative alleles
        return ",".join(sorted(alternatives))  # Sorted for consistent ordering
    
    # Apply this function row-wise
    merge_vcf_imp_exi["ALT"] = merge_vcf_imp_exi.apply(lambda row: fill_alt_from_y(row, samp_y_cols), axis=1)

    # print(samp_list)
    # print(merge_vcf_imp_exi)
    #COMBINE with new SNPs
    ref_dic = ref_dic[["Ref",args.vcf_ref,"Anc","Der"]]
    ref_dic = ref_dic.drop_duplicates()
    # print(ref_dic)
    #merge
    merge_vcf_imp_exi = merge_vcf_imp_exi.merge(ref_dic, left_on="POS", right_on=str(args.vcf_ref), how= "left")
    # print(merge_vcf_imp_exi)
    def normalize(x):
        return str(x).upper()

    mask = merge_vcf_imp_exi["#CHROM"] == args.vcf_chr

    for idx in merge_vcf_imp_exi[mask].index:
        ref = normalize(merge_vcf_imp_exi.at[idx, "REF"])

        # Preserve original ALT order
        alt_list = [
            normalize(a)
            for a in str(merge_vcf_imp_exi.at[idx, "ALT"]).split(",")
            if a
        ]

        alt_set = set(alt_list)  # for fast membership checking

        for col in ["Anc", "Der"]:
            val = normalize(merge_vcf_imp_exi.at[idx, col])

            if val != ref and val not in alt_set:
                alt_list.append(val)
                alt_set.add(val)

        merge_vcf_imp_exi.at[idx, "ALT"] = ",".join(alt_list)

    merge_vcf_imp_exi = merge_vcf_imp_exi.iloc[:, :-4]

    # print(merge_vcf_imp_exi)
    
    # print(ee)
    
    for samp in samp_list:
        # print(samp)
        # Apply the function row-wise
        merge_vcf_imp_exi[f"{samp}_y"] = merge_vcf_imp_exi.apply(lambda row: get_variant_code(row, samp), axis=1)
        # print(merge_vcf_imp_exi)
        # print(ee)
        
        merge_vcf_imp_exi.loc[merge_vcf_imp_exi[f"{samp}_x"].astype(str).str.startswith(".") & merge_vcf_imp_exi[f"{samp}_y"].notna(), f"{samp}_x"] = merge_vcf_imp_exi[f"{samp}_y"]
        
        # merge_vcf_imp_exi.loc[merge_vcf_imp_exi[f"{samp}"].astype(str).str.startswith(".") & merge_vcf_imp_exi[f"{samp}_y"].notna(), f"{samp}_x"] = merge_vcf_imp_exi[f"{samp}_y"]
        # print("E")
        # print(merge_vcf_imp_exi)
        # print(ee)
        merge_vcf_imp_exi = merge_vcf_imp_exi.drop(columns=[f"{samp}_y"])
        merge_vcf_imp_exi = merge_vcf_imp_exi.rename(columns={f"{samp}_x":f"{samp}"})
    
    old_col = f"{args.vcf_ref}_x"
    new_col = args.vcf_ref

    if old_col in merge_vcf_imp_exi.columns:
        merge_vcf_imp_exi = merge_vcf_imp_exi.rename(columns={old_col: new_col})
    # print(merge_vcf_imp_exi)
    # print(ee)
    
    
    
    merge_vcf_imp_exi = merge_vcf_imp_exi.drop(columns=[args.vcf_ref,"marker"])
    merge_vcf_imp_exi = merge_vcf_imp_exi.drop_duplicates(subset="POS").reset_index(drop=True) 
    
    df_dic_merge_new_ref = df_dic_merge_new.merge(ref_dic, left_on=args.vcf_ref, right_on=args.vcf_ref, how="left")
    
    #Add to Alt column the observations in the imputed dataset
    df_dic_merge_new_ref = df_dic_merge_new_ref.drop_duplicates(subset = [args.vcf_ref])
    df_dic_merge_new_ref["Alt"] = df_dic_merge_new_ref.apply(
    lambda row: ",".join(
        sorted(set(
            val for val in row.iloc[2:-2]
            if pd.notna(val) and val != row["Ref"]
        ))
    ) if any(pd.notna(val) and val != row["Ref"] for val in row.iloc[2:-2]) else ".",
    axis=1
    )
    
    # print(df_dic_merge_new_ref)
    #replace N with .:.
    for samp in samp_list:
        df_dic_merge_new_ref[samp] = df_dic_merge_new_ref[samp].replace("N", ".:.").fillna(".:.")    
        #add imputation code
        df_dic_merge_new_ref[f"{samp}"] = df_dic_merge_new_ref.apply(
            lambda row: (
                "0:999" if row[f"{samp}"] == row["Ref"]
                else row[f"{samp}"] if row[f"{samp}"] == ".:."
                else (
                    # Attempt to find match in ALT list
                    (lambda alt_alleles, val: (
                        f"{alt_alleles.index(val) + 1}:999" if val in alt_alleles else val
                    ))(str(row["Alt"]).split(","), row[f"{samp}"])
                ) if pd.notna(row[f"{samp}"])
                else np.nan
            ),
            axis=1
        )

    # print(df_dic_merge_new_ref)
    df_dic_merge_new_ref = df_dic_merge_new_ref.drop(columns=["marker"])
    df_dic_merge_new_ref = df_dic_merge_new_ref.rename(columns={args.vcf_ref:"POS"})
    df_dic_merge_new_ref.insert(0, "#CHROM", args.vcf_chr)
    
    df_dic_merge_new_ref.insert(2, "FORMAT", ".")
    df_dic_merge_new_ref.insert(2, "INFO", ".")
    df_dic_merge_new_ref.insert(2, "FILTER", ".")
    df_dic_merge_new_ref.insert(2, "QUAL", ".")
    
    alt_col = df_dic_merge_new_ref.pop("Alt")
    df_dic_merge_new_ref.insert(2, "ALT", alt_col)
    
    ref_col = df_dic_merge_new_ref.pop("Ref")
    df_dic_merge_new_ref.insert(2, "REF", ref_col)
    df_dic_merge_new_ref.insert(2, "ID", ".")
    
    # print(df_dic_merge_new_ref)
    # print(merge_vcf_imp_exi)
    
    #combine dfs
    vcf_df = pd.concat([merge_vcf_imp_exi, df_dic_merge_new_ref], ignore_index=True)
    vcf_df ["POS"] = vcf_df ["POS"].astype(int)
    vcf_df  = vcf_df.sort_values(by="POS")
    vcf_df = vcf_df.drop(columns=["Anc","Der"])

    
    #Now replace the imputed content in the multisample vcf file
    
    header_lines = []
    output_path2 = os.path.join(args.output + "_phyloimputed.vcf")

    # Read gzip VCF header (not the data)
    with gzip.open(output_path, "rt") as infile:
        for line in infile:
            if line.startswith("#CHROM"):
                header_lines.append('##FORMAT=<ID=IMP,Number=1,Type=String,Description="Imputation flag: IMP=999 means phyloimputed">\n')
                break
            else:
                header_lines.append(line)

    # Write header lines + DataFrame header + data rows to uncompressed VCF file
    with open(output_path2, "wt") as outfile:
        # Write the original header lines
        for line in header_lines:
            outfile.write(line)

        # Write the DataFrame header line (tab-separated column names)
        outfile.write('\t'.join(vcf_df.columns) + '\n')

        # Write DataFrame rows
        for _, row in vcf_df.iterrows():
            outfile.write('\t'.join(map(str, row)) + '\n')
        
    if os.path.exists(output_path):
        os.remove(output_path)
    with open(log_filename, "a") as f:
        f.write(f'[{timestamp}]\tFinished saving phyloimputed multi-sample vcf file in {output_path}\n')

    
if args.nucleotide:
    with open(log_filename, "a") as f:
        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        f.write(f'[{timestamp}]\tStart converting allelic states to nucleotides\n')
    dic_2 = dic.drop(columns=["GRCh37","GRCh38","T2T","Hg"])
    # dic_2 = dic_2.drop_duplicates()
    dic_2 = dic_2.drop_duplicates(subset='marker')
    dic_indexed = dic_2.drop_duplicates(subset='marker').set_index('marker')

    # first define the transformation function
    def transform_row(row):
        marker = row.name
        if marker not in dic_indexed.index:
            return None
        anc = dic_indexed.at[marker, 'Anc']
        der = dic_indexed.at[marker, 'Der']
        def convert(val):
            if val == 'a':
                return anc.lower()
            elif val == 'A':
                return anc.upper()
            elif val == 'd':
                return der.lower()
            elif val == 'D':
                return der.upper()
            elif val == 'X':
                return 'N'
            else:
                return "unexpected entry"  # leave unchanged if unexpected
        return row.apply(convert)
    #apply function
    df = df.apply(transform_row, axis=1)
    df = df.dropna()  # Drop rows that returned None
    path_output = args.output + "_phyloimputed_nucl.csv"
else:
    path_output = args.output + "_phyloimputed.csv"

#Save outputs
path_conflicting_SNPs = args.output + "_conflicting_SNPs.csv"
haplogroup_info = args.output + "_haplogroups.csv"

df.to_csv(path_output, index=True, sep='\t')
np.savetxt(haplogroup_info, Haplogroup_info, delimiter="\t", fmt="%s", comments="")
np.savetxt(path_conflicting_SNPs, questionable_SNPs, delimiter="\t", fmt="%s", comments="")


end_time = time.time()
elapsed_time = end_time - start_time
with open(log_filename, "a") as f:
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    f.write(f'[{timestamp}]\tSaved output files\n')
    f.write(f'[{timestamp}]\tFinished run\n')
    f.write(f"""
                  ┌─ PHYLOIMPUTE COMPLETE!
              ┌───┤
              │   └─ Script completed in {elapsed_time:.2f} seconds.
──────────────┤
              └─ :)
""")


print("\n\n")

print(f"""
                  ┌─ PHYLOIMPUTE COMPLETE!
              ┌───┤
              │   └─ Script completed in {elapsed_time:.2f} seconds.
──────────────┤
              └─ :)
""")
print("\n")

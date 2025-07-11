#don't use commas in snp names
import argparse
import time
import shutil
import glob
import os, sys
import pandas as pd
import math
import numpy as np
import subprocess
import warnings

#User input paths
parser = argparse.ArgumentParser()

parser.add_argument("-input_format", choices=["vcf","csv"],help="Option 'vcf' is used when the available input format are .vcf or vcf.gz files. The vcf files need to be quality filtered already! The 'csv' option is used for tab-delimited csv file of a dataframe with individuals as columns and variants as rows. "
  "Insert variant data as 'A' for ancestral, 'D' for derived allele, and 'X' for missing data. Add variant names as first column and sequence name as header row. "
  "For both options: Refrain from using comma separated marker or sequence names.")

parser.add_argument("-input", required=True, help="For 'vcf' give the folder to the vcf files. For 'csv', give the path to one csv file. If -freqmap is specified, this needs to specify phyloimputed (or external) tab-separated csv file for allele frequency map.")
parser.add_argument("-output", required=True, help="Path to output folder.")
parser.add_argument("-vcf_ref", choices=["GRCh37","GRCh38","T2T"],help="When using the 'vcf' option, give the version of the reference genome used for the alignment: GRCh37, GRCh38, or T2T") 
parser.add_argument("-vcf_chr", help="How is the chromosome you want to analyze, named in your vcf file? This depends on the reference genome used and can be e.g. 'Y','chrY', 'NC_000024.9' for the human Y chromosome.")
parser.add_argument("-vcf_dic", help="If you are using a custom tree, you need to provide the path to a tab-separated .csv dictionary file that gives information on the genetic markers of interest that are also provided in the custom tree. Using the following column names: 'marker'(=marker names that are identical with the custom tree marker names),'GRCh37' and/or 'GRCh38' and/or 'T2T','Anc'(=ancestral allele), 'Der'(=derived allele)")
parser.add_argument("-tree", choices=["Y_minimal","NAMQY","ISOGG_2020"],help="Optional: path to tab-separated custom file of the phylogenetic SNP tree.")
parser.add_argument("-customtree", help="The path to a custom tree has to be provided in the same format as the tree files Y_minimal.csv, etc. Custom tree need to start with a root snp, e.g. 'ROOT'.")

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

#GENERATE ALLELE FREQUENCY MAP


if args.freqmap:
    # if type(args.color) != str: #default
        # args.color = "blue"
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
        # cmd.append('-country', args.country)  
        cmd.extend(['-country'] + args.country) 
    # print(cmd)
    subprocess.run(cmd)    

    sys.exit()


if not args.input_format:
    print("Please define '-input_format' to run PhyloImpute.")
    sys.exit()
    
#FOR VCF INPUT FILE
print("PhyloImpute v1.2 is running...")
def check_for_vcf_files(directory):
    vcf_files = glob.glob(os.path.join(directory, '*.vcf'))
    if vcf_files:
        print(f'\nFound .vcf files:')
        print(vcf_files)
        return vcf_files
    else:
        print('\nNo .vcf files found.')
        return []

#get_marker_info_from_vcf by using a dictionary file
if args.input_format == "vcf":
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
        print("isogg")
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
    #Make input table for Phyloimpute
    #turn positions to strings, so they can be recognized later in the data tables
    input_table = pd.DataFrame({'pos':dic_pos, 'Sample':dic_marker})
    column_missing_data = len(input_table)*["X"]


def extract_sample_name(file_path):
    split_string = file_path.split("/")
    sample_name_format = split_string[-1]
    sample_name = sample_name_format.replace(".csv", "")
    sample_name = sample_name_format.replace(".vcf", "")
    # print(sample_name)
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
    #Filter chromosome column
    df = df[df["#CHROM"].isin([args.vcf_chr])]
    #Filter positions of variants from database
    df = df[df["POS"].isin(dic_pos)]
    #Also filter that is the derived allele that we see:
    ##Make list of variant positions in the same order
    vcf_pos_list = df["POS"].tolist()
    dic_vcf_pos = dic[dic[args.vcf_ref].isin(vcf_pos_list)]


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
    df["der_from_dic_for_pos"] = sorted_dic_vcf_pos["Der"]
    #compare two columns in a dataframe, if they dont match, remove the row
    df_filtered = df[df['ALT'] == df['der_from_dic_for_pos']]
    #Make list of the position column
    df_filtered_list = df_filtered["POS"].values.tolist()
    return df_filtered_list



if args.input_format == "vcf":
    vcf_files = check_for_vcf_files(args.input)
    if type(args.vcf_ref)==str and type(args.vcf_chr)==str:
        print("\nNecessary input files are provided!\n")
        #Now that checking is over, we loop through 
        for vcf_file in vcf_files:
            print(f'Processing file: {vcf_file}')
            derived_variant_pos=filter_vcf_file(vcf_file)
            # print(derived_variant_pos)
            sample_name = extract_sample_name(vcf_file)
            #Add sample column to input table with "X"s in it
            input_table[sample_name] = column_missing_data
            #filter for the positions that the sample has derived alleles in the vcf file
            input_table_D = input_table[input_table["pos"].isin(derived_variant_pos)]
            input_table_D.loc[:,sample_name]="D"
            input_table = pd.concat([input_table_D, input_table[~input_table["pos"].isin(derived_variant_pos)]])
        #replace pos column and reindex
        input_table.reset_index(drop=True, inplace=True)
        input_table["pos"]=range(1,len(input_table)+1)
        df = input_table
        df.set_index("Sample", inplace=True)
        #rename column pos
        df.rename(columns={'pos':0}, inplace=True)
    else:
        print("\nOne or both of the input files vcf_ref or vcf_chr are missing! Please provide and rerun. See help function for more information. \n")



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

if type(args.tree) == str:
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

Haplogroup_info = ["Sample\tHg\tConfidence_value\tPenalty_value1\tPenalty_value2"]
for col_name, col in df.items(): #iteritems()
    lengths_list = []
    filtered_D= df[col_name]=="D"
    sample_D_list = df[filtered_D].index.tolist()#make a list of the index values (SNP names) of the filtered dataframe
    filtered_A= df[col_name]=="A"
    sample_A_list = df[filtered_A].index.tolist()#make a list of the index values (SNP names) of the filtered dataframe
    sample_A_set=set(sample_A_list)
    #Compare each sample's derived lists with the database lists
    for db_list in tree_lists:
        processed_D_list=preprocess_set(sample_D_list)
        processed_db_list=preprocess_set(db_list)
        overlapping_entries=processed_D_list.intersection(processed_db_list)
        len_overlap = len(overlapping_entries)
        lengths_list.append(len_overlap)
    max_index = lengths_list.index(max(lengths_list))
       
    #get the database row that will be used to replace missing data
    #are any of these items overlapping with ancestral observations in the sample?
    new_D_db =set(tree_lists[max_index])
    processed_set1=preprocess_set(sample_A_list)
    max_d=preprocess_set(tree_lists[max_index])
    
    #GET INFO ON HG OF SAMPLE           ##HG
    overlapping_derived=processed_D_list.intersection(max_d) #overlap between derived snps and snps of max branches
    
    if len(overlapping_derived) == 0:
        Hg_sample = "No hg pred possible, because no informative derived alleles"
        Hg_support = "-"
        Hg_penalty1 = "-"
        Hg_penalty2 = "-"
    else:
        Hg_sample = hg_table[max_index]
        Hg_support = str(len(overlapping_derived))+"/"+str(len(max_d))
        #Penalty comprises of variants in the max_d_branch that are ancestral (but are expected to be derived)
        #But also of variants that are derived but not in the max_d_branch
        Hg_penalty_part1 = processed_set1.intersection(max_d)#How many are ancestral in our max res branch
        Hg_penalty1 = str(len(Hg_penalty_part1))+"/"+str(len(max_d))#How many are ancestral
        
        Hg_penalty_part2 = processed_D_list - max_d #How many in other branches are derived #TO DO: exclude SNPs that are not in tree
        Hg_penalty_part2 = Hg_penalty_part2 - df_exclusive_SNPs #remove SNPs that are derived in data, but not in tree
        if "ROOT" in tree_snps:
            total_snps_tree = len(tree_snps)-2 #Remove title and ROOT
        else:
            total_snps_tree = len(tree_snps)-1 #Remove title, No ROOT
        Hg_penalty2 = str(len(Hg_penalty_part2))+"/"+str(total_snps_tree) #divided by total number of variants
        
    Haplogroup_info.append(col_name+"\t"+Hg_sample+"\t"+str(Hg_support)+"\t"+str(Hg_penalty1)+"\t"+str(Hg_penalty2)) #col_name is the sample name

    
    overlapping_entries=processed_set1.intersection(max_d)
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
            # print(parallel_markers)
            #remove nan
            parallel_markers = {x for x in parallel_markers if pd.notnull(x)}
            parallel_markers = {x for x in parallel_markers if x is not np.nan}
            # print(parallel_markers)
            # print(ee)
            A_list.append(list(parallel_markers))
            #to set, to remove duplicate entries   
    A_list=set([item for sublist in A_list for item in sublist]) #flatten list in lists
    A_list = A_list - set(tree_lists[max_index]) #set(tree_lists[max_index]) are the derived snps
    # print(A_list[:-6])
    A_list_markers=preprocess_set(A_list)
    #also remove the derived single (already processed snps) that are derived in the sample:
    A_list_markers = A_list_markers -set(sample_D_list)
    #now remove the ones that are already known to be ancestral
    A_list_markers=A_list_markers-sample_A_set

    if len(A_list_markers) > 0:
        filtered_df = df[df.index.isin(A_list_markers)]
        unfiltered_df = df[~df.index.isin(A_list_markers)] #ALL OTHER VARIANTS to combine in the end with the altered filtered_df 
        filtered_df.loc[:,col_name]="a" #make the all derived (filling gaps step) #*for inference
        df=pd.concat([filtered_df,unfiltered_df])

#flatten list of sets:
questionable_SNPs=list(questionable_SNPs)

#Save outputs
path_output = args.output + "_phyloimputed.csv"
path_conflicting_SNPs = args.output + "_conflicting_SNPs.csv"
haplogroup_info = args.output + "_haplogroups.csv"

df.to_csv(path_output, index=True, sep='\t')
np.savetxt(haplogroup_info, Haplogroup_info, delimiter="\t", fmt="%s", comments="")
np.savetxt(path_conflicting_SNPs, questionable_SNPs, delimiter="\t", fmt="%s", comments="")


print("\nPHYLOIMPUTE COMPLETE!")

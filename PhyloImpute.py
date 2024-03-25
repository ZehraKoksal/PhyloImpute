#!/mnt/ngs/scratch_areas/nxd426/python3.6_env/bin/python

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

#User input paths


parser = argparse.ArgumentParser()
parser.add_argument("input", help="Tab-delimited csv file of a dataframe with individuals as columns and variants as rows. "
  "Insert variant data as 'A' for ancestral, 'D' for derived allele, and 'X' for missing data. Add variant names as first column and sequence name as header row. "
  "Refrain from using comma separated marker or sequence names.")
parser.add_argument("output", help="Path to output folder.")
parser.add_argument("--tree", choices=["Y_minimal","Y_ISOGG_2019-20","Y_NAM-Q-Y"],help="Optional: path to tab-separated custom file of the phylogenetic SNP tree.")
parser.add_argument("--customtree", help="The path to a custom tree has to be provided in the same format as the tree files Y_minimal.csv, etc. Custom tree need to start with a root snp, e.g. 'ROOT'.")

args = parser.parse_args()

current_dir = os.path.dirname(os.path.abspath(__file__))

#import sample dataframe and variant order 
df = pd.read_csv(args.input, sep="\t", index_col=1)
# df = pd.read_csv("/mnt/ngs/scratch_areas/nxd426/3_Q_AmpliSeq/Yleaf/Final_Yleaf_RD20_new_names/OVERVIEW_Qpanel_422_samples.csv", sep="\t", index_col=1)
#skip numbering column
df = df.iloc[:,1:]
#replace missing data with X and (*) with X
df = df.fillna("X")
df = df.replace(to_replace=r'\(.*\)', value='X',regex=True)


#TREE to DATABASE
# if args.tree == "custom":
if type(args.customtree) == str:
    print("E")
    print(args.customtree)
    # custom_tree_path = os.path.join(current_dir, "Database_generator.py")
    tree = pd.read_csv(args.customtree, sep="\t", header=None)
    print("TREE")


    for col_name, col in tree.iteritems():
        first_non_nan_index = col.first_valid_index() ##Get the index with the first non nan
        for index, value in col.iloc[first_non_nan_index:].items():#tree.iterrows():
            if index==len(tree)-1:
                break
            if pd.notna(index):
                marker=value
            if col_name==0:
                tree.loc[:,col_name]  = col.fillna(marker) #replace all nan with the first marker non-nan
                break #so we can skip to the next column
            elif pd.isna(col.iloc[index+1]): #and index<len(tree)-2:#pd.isna(index+1,col):
                prev_upstream=tree.loc[index,col_name-1]
                if tree.loc[index+1,col_name-1]==prev_upstream :
                    tree.loc[index+1,col_name]  = marker
                    if index==len(tree)-2:
                        break
    custom_tree_path = current_dir + "/custom.csv"
    tree.to_csv(custom_tree_path, index=False, sep='\t')
    
    # tree.to_csv("/mnt/ngs/scratch_areas/nxd426/Missing_data_phylogeny/Q_panel_tree_SNP_lib.csv", index=False, sep='\t')

    print("SNP database generated and saved!")
    # if os.path.exists(custom_tree_path):
        # subprocess.run(["python",custom_tree_path]) #open the database generator python


if type(args.tree) == str:
    filename = args.tree + ".csv"
    tree_csv_path = os.path.join(current_dir,filename)
    # tree = pd.read_csv("/mnt/ngs/scratch_areas/nxd426/Missing_data_phylogeny/Q_panel_tree_SNP_lib.csv", sep="\t")#, index_col=1)
    tree = pd.read_csv(tree_csv_path, sep="\t")


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

# Add all snp names that are not in the sample df, as rows filled with X's
df_snps = set(df.index.tolist())

adding_snps = tree_snps - df_snps
adding_snps = {x for x in adding_snps if x is not np.nan}
adding_snps = {x for x in adding_snps if x != ""}

#add the missing snps from the tree to the data
df_adding_snps = pd.DataFrame(index=adding_snps, columns=df.columns, data="X")
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
for col_name, col in df.iteritems():
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
    #get the database row that will be used to replace mising data

    #are any of these items overlapping with ancestral observations in the sample?
    new_D_db =set(tree_lists[max_index])
    processed_set1=preprocess_set(sample_A_list)
    max_d=preprocess_set(tree_lists[max_index])
    overlapping_entries=processed_set1.intersection(max_d)
    if len(overlapping_entries)>0:
        questionable_entry = str(col_name)+":"+str(overlapping_entries)
        questionable_SNPs.append(questionable_entry)

    
    #A)get the SNPs from the database row, but exclude thoes that were ancestral in the sample
    to_replace_derived = max_d.difference(processed_set1)
    new_D_db = max_d-processed_set1-processed_D_list #the snps in the database row excluding the ancestral and derived variants from the current sample
    if len(new_D_db) > 0:
        filtered_df = df[df.index.isin(new_D_db)]
        unfiltered_df = df[~df.index.isin(new_D_db)] #ALL OTHER VARIANTS to combine in the end with the altered filtered_df 
    # filtered_df = df[df.index.map(lambda x: any(substring in x.split(',') or substring == x for substring in to_replace_derived))]
        filtered_df[col_name]="d" #make the all derived (filling gaps step) #*for inference
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

    if len(A_list_markers) > 0:
        filtered_df = df[df.index.isin(A_list_markers)]
        unfiltered_df = df[~df.index.isin(A_list_markers)] #ALL OTHER VARIANTS to combine in the end with the altered filtered_df 
        # print(filtered_df)
        # print(unfiltered_df)
        filtered_df[col_name]="a" #make the all derived (filling gaps step) #*for inference
        df=pd.concat([filtered_df,unfiltered_df])

#flatten list of sets:
# questionable_SNPs = {item for s in questionable_SNPs for item in s}
questionable_SNPs=list(questionable_SNPs)

df.to_csv("/mnt/ngs/scratch_areas/nxd426/Missing_data_phylogeny/output_imputed_data.csv", index=True, sep='\t')
np.savetxt("/mnt/ngs/scratch_areas/nxd426/Missing_data_phylogeny/questionable_SNPs.csv", questionable_SNPs, delimiter="\t", fmt="%s", comments="")

print("PHYMPUTE COMPLETE!")
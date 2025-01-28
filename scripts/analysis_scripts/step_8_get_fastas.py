#this script will take a csv file and make three fasta files needed to align and remove redundancy
#1. pymol.fasta- the fasta sequence from the pymol structure: the index here will match the reported index in the csv file
#2. pdb.fasta- the fasta sequences pulled from pdb
#3. reference.fasta- non-redundant reference sequences from the pdb.fasta

##I am writing this script to be partnered with 2_adjust_index.py but I have pymol installed on my computer and needle on the bioinfo conda env so I can't have them in the same script

import argparse
import os
import pandas as pd
import numpy as np
from pymol import cmd, stored
import requests
import re

parser = argparse.ArgumentParser(description='This script will make three fasta files needed to align and remove redundancy')
parser.add_argument('-i', '--inputcsv', help='non-redundant csv file', required=True)
parser.add_argument('-ft','--fasta_type', help='type of fasta file to make, i.e. (pymol) pymol, (pdb) pdb, or (ref) reference', required=True)
parser.add_argument('-pdb','--pdbfasta', help='fasta file from reference sequence, just one', required=False)
parser.add_argument('-o','--output', help='output fasta file', required=True)
args= parser.parse_args()

#read in csv file
df = pd.read_csv(args.inputcsv)

#make a fasta file from the pymol structure in two functions
def get_pymol_seq(pdb_id, chain_name):
    #fetch pdb file by id
    cmd.fetch(pdb_id)
    #print("step1")
    #if cmd.fetch fails, return "pymol fail"
    if cmd.count_atoms() == 0:
        return "pymol fail"
    #iterate through chain and save residue number and residue name as a list
    stored.residues = []
    cmd.iterate('chain %s' % chain_name, 'stored.residues.append([resi, resn])')
    cmd.delete(pdb_id)
    #delete stored .cif file
    os.remove(pdb_id+".cif")
    print(f"pymol progress {pdb_id}_{chain_name}")
    #make a new list where each residue is only listed once
    fasta_wmg = []
    for i in stored.residues:
        if i not in fasta_wmg:
            fasta_wmg.append(i)
    #print("step2")
    #make a new list where all residues are digits
    fasta_digit = []
    for i in fasta_wmg:
        if i[0].isdigit() == True:
            fasta_digit.append(i)
    #print("step3")
    #make a new list where only A C G U is included
    fasta = []
    for i in fasta_digit:
        #if item[0] is not a number, remove i from fasta list
        if i[1]== "A" or i[1] == "C" or i[1] == "G" or i[1] == "U":
            fasta.append(i)
    #print("step4")
    #if there i[0] = 1 is anywhere in fasta, remove all indexes before it
    first_nt = None
    for index, i in enumerate(fasta):
        if i[0] == '1':
            first_nt = index
    if first_nt != None:
        fasta = fasta[first_nt:]
    #print("step5")
    #if fasta i[0] doesn't start with one, add ['1','N'] to the beginning of the list
    if fasta[0][0] != '1':
        fasta.insert(0,['1','N'])
    #print("step6")
    #order residues numerically by index number
    fasta = sorted(fasta, key=lambda x: int(x[0]))
    #remove any item with an index <1
    for index, item in enumerate(fasta):
        if int(item[0]) < 1:
            fasta.pop(index)
    #print("step7")
    #insert new list items to make sure the residues are sequential
    for index, item in enumerate(fasta):
        #check if index number is equal to list item
        if int(item[0]) == index+1:
            continue
        #if it isn't, insert list item at index number
        elif int(item[0]) != index+1:
            fasta.insert(index,[str(index+1), "N"])
    #print("step8")
    #reformat for a fasta file
    seq_string = ""
    for i in fasta:
        seq_string += i[1]
    fasta = seq_string
    return(fasta)


def pymol_to_fasta(res_info_df):
    #save unique data_ext values to pdb_chain_list
    pdb_chain_list = res_info_df["data_ext"].unique()
    #split pdb_chain_list into a list of lists with pdb_id and chain_id using pandas
    pdb_chain_list = [i.split("_") for i in pdb_chain_list]
    #iterate through pdb_chain_list and get the pymol_seq for each pdb_id and chain_id, save to dataframe
    fasta_list = []
    count = 0
    for i in pdb_chain_list:
        pdb_id = i[0]
        chain_id = i[1]
        pdb_chain = ">"+pdb_id + "_" + chain_id
        fasta = get_pymol_seq(pdb_id, chain_id)
        fasta_list.append([pdb_chain, fasta])
        print(f"progress {count+1}/{len(pdb_chain_list)}")
        count += 1
    #write fasta_list to a file
    with open(str(args.output) + ".fasta", "w") as f:
        for i in fasta_list:
            f.write(i[0] + "\n" + i[1] + "\n")


def get_pdb_seq(pdb_id, chain_name):
    print(f"pdb progress {pdb_id}_{chain_name}")
    #define url
    url = "https://www.rcsb.org/fasta/entry/" + pdb_id + "/download"
    response = requests.get(url)
    if response.status_code == 200:
        full_fasta = response.text
    else:
        return (f"error with {pdb_id}_{chain_name} download")
    #save each line as an item to a list
    full_fasta = full_fasta.splitlines()
    fasta_list = []
    for line in full_fasta:
        if line.startswith(">"):
            header = line
            fasta_list.append([header, ""])
        else:
            fasta_list[-1][1] += line
    #filter the fasta to remove proteins
    filtered_fasta = [item for item in fasta_list if "V" not in item[1]]
    filtered_fasta2 = [item for item in filtered_fasta if "L" not in item[1]]
    #keep only the list item that contains "auth chain_name" but doesn't contain "protein" in the first index
    auth_label = "[auth " + chain_name +"]"
    chain_label = "|Chain " + chain_name + "|"
    space_chain = " " + chain_name
    fasta_out = []
    for i in filtered_fasta2:
        if auth_label in i[0]:
            fasta_out = i[1]
            break
        elif chain_label in i[0]:
            fasta_out = i[1]
            break
        elif re.search(r'Chains\b[^|]*' + re.escape(space_chain), i[0]):
            fasta_out = i[1]
            break
        else:
            fasta_out= [">" + pdb_id + "_" + chain_name, "pdb fail"]
    return fasta_out


def pdb_to_fasta(res_info_df):
    #save unique data_ext values to pdb_chain_list
    pdb_chain_list = res_info_df["data_ext"].unique()
    #split pdb_chain_list into a list of lists with pdb_id and chain_id using pandas
    pdb_chain_list = [i.split("_") for i in pdb_chain_list]
    #iterate through pdb_chain_list and get the pdb_seq for each pdb_id and chain_id, save to dataframe
    fasta_list = []
    count = 0
    for i in pdb_chain_list:
        pdb_id = i[0]
        chain_id = i[1]
        pdb_chain = ">"+pdb_id + "_" + chain_id
        fasta = get_pdb_seq(pdb_id, chain_id)
        #if the length of the fasta is equal to "Chain_length_reference", add to fasta_list
        if len(fasta) == res_info_df[res_info_df["data_ext"]==pdb_id + "_" + chain_id]["Chain_length_reference"].values[0]:
            fasta_list.append([pdb_chain,fasta])
        else:
            fasta_list.append([pdb_chain,"length fail"])
        print(f"progress {count+1}/{len(pdb_chain_list)}")
        count += 1
    #write fasta_list to a file
    with open(str(args.output) + ".fasta", "w") as f:
        for i in fasta_list:
            f.write(i[0] + "\n" + i[1] + "\n")


def get_ref_seq(res_info, pdb_fast):
    #add a new column to the res_info dataframe with rna_org
    res_info["rna_org"] = res_info["Molecule"] + "_" + res_info["Source_Organism"]
    #make a list of unique rna_org values
    rna_org_list = res_info["rna_org"].unique()
    #for each value of rna_org, get the average Chain_length_reference where value == rna_org in res_info
    ref_length = []
    for i in rna_org_list:
        mean = res_info[res_info["rna_org"]==i]["Chain_length_reference"].mean()
        ref_length.append([i,mean])
    #get the data_ext value and rna_org where Chain_length_reference is the longest but not more than 1.5 times the mean
    ref_list = []
    for i in ref_length:
        max_length = res_info[(res_info["rna_org"]==i[0]) & (res_info["Chain_length_reference"]<i[1]*1.5)]["data_ext"].max()
        #append max_length and i[0] to ref_list
        ref_list.append([i[0], max_length])
    #read in pdb_fast as a dictionary where the line starting with > is the key and the following line is the value
    with open(pdb_fast, "r") as f:
        pdb_dict = {}
        for line in f:
            if line.startswith(">"):
                key = line.strip()
            else:
                pdb_dict[key] = line.strip()
    #make a data frame with the ref_list and pdb_dict value where ref_list == key
    ref_list2 = []
    for i in ref_list:
        pdb_id = ">" + str(i[1])
        name = ">" + str(i[0]) + "|" + str(i[1])
        ref_list2.append([name,pdb_dict[pdb_id]])
    #save ref_list2 to a file
    with open(str(args.output) + ".fasta", "w") as f:
        for i in ref_list2:
            f.write(i[0] + "\n" + i[1] + "\n")





#main function
#pymol_to_fasta(df)




if args.fasta_type == "pymol":
    pymol_to_fasta(df)
if args.fasta_type == "pdb":
    pdb_to_fasta(df)
if args.fasta_type == "ref":
    #if args.pdbfasta isn't provided, print an error
    if args.pdbfasta == None:
        print("please provide a pdb fasta file")
    else:
        get_ref_seq(df, args.pdbfasta)

#this script takes the 3 fasta files from 1_get_fastas.py as well as the original csv 
#it adjust the residue index to one common fasta (reference.fasta) within species for rRNAs
#it's written to easily be combined with 1_get_fastas.py but since they are installed weird on my computer I kept them separate

#the general flow of the steps are:
#1. align pymol fasta to pdb fasta to fill in the gaps (this can be later modified downstream to get flanking sequences) and get pdb adjusted index
#2. align pdb fasta to reference fasta to get reference adjusted index
#3. output a .csv file with pdb, chain, res1 index, res2 index, ref pdb_chain, res1 ref adjusted index and res2 ref ajusted index

import argparse
import pandas as pd
import csv
import os
from Bio.Emboss.Applications import NeedleCommandline
from Bio.Application import ApplicationError
import subprocess
#from Bio.Align.Applications import NeedleCommandline


#set arguments
parser = argparse.ArgumentParser(description='This script will make three fasta files needed to align and remove redundancy')
parser.add_argument('-i', '--inputcsv', help='non-redundant csv file', required=True)
parser.add_argument('-py','--pymolfasta', help='fasta file from pymol', required=True)
parser.add_argument('-pdb','--pdbfasta', help='fasta file from pdb', required=True)
parser.add_argument('-ref','--referencefasta', help='fasta file from reference', required=True)
parser.add_argument('-o','--output', help='output fasta file', required=True)
args= parser.parse_args()

#read in csv file
df = pd.read_csv(args.inputcsv)
#add a new column "pdb_chain" that takes the first item when data_ext is split by "."
df["pdb_chain"] = df["data_ext"].str.split(".").str[0]
#store the molecule and source organism chain as ref_org
df["mol_org"] = ">" + df["Molecule"] + "_" + df["Source_Organism_chain"]



#define function to read in fasta file as dictionary
def fasta_to_dict(fasta_file):
    fasta_dict = {}
    header = None
    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line
                fasta_dict[header] = ""
            elif header != None:
                fasta_dict[header] += line
    return fasta_dict

pymol_dict = fasta_to_dict(args.pymolfasta)
pdb_dict = fasta_to_dict(args.pdbfasta)
reference_dict = fasta_to_dict(args.referencefasta)



#define function to run needle alignment and reformat output to be useable
#def run_alignment(fasta1,fasta2):

#run needle command
def run_needle(asequence, bsequence):
    #if any character in asequence is not A, C, G, U or N, replace it with N
    asequence = "".join([i if i in ["A", "C", "G", "U", "N"] else "N" for i in asequence])
    bsequence = "".join([i if i in ["A", "C", "G", "U", "N"] else "N" for i in bsequence])

    try:
        needle_cline = NeedleCommandline(
            asequence="asis::" + asequence,
            bsequence="asis::" + bsequence,
            gapopen=10,
            gapextend=0.5,
            outfile="align.txt")
        stdout, stderr = needle_cline()


        with open("align.txt", "r") as alignment_file:
            alignment_text = alignment_file.read()
        return alignment_text
        
    except ApplicationError as e:
        #write "alignment failed: " + str(e) to align.txt
        with open("align.txt", "w") as f:
            f.write("alignment failed: " + str(e))


#print(run_needle("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUCCUCGUNNCCCAANGGNCACGGCNNCUGGCUNCGAACCAGAAGANUNCAGGNNCANGUCCUGGCGGGGAAGCCA","UUCCUCGUGGCCCAAUGGUCACGGCGUCUGGCUNCGAACCAGAAGAUUCCAGGUUCAAGUCCUGGCGGGGAAGCCA"))



#reformat alignment output
def alignment_format(alignment_file):
    #read in alignment output
    with open(alignment_file, "r") as f:
        alignment_text = f.read()
    #if alignment_text starts with "alignment failed", return "alignment failed"
    if alignment_text.startswith("alignment failed"):
        return "alignment failed"
    align_lines = []
    asequence = []
    bsequence = []
    #split alignment output into lines
    for line in alignment_text.splitlines():
        #if line starts with "asis"
        if line.startswith("asis"):
            #split line by space
            line = line.split()
            #concatonate 3rd element of line to asequence
            align_lines.append(line[2])
    #store lines in appropriate list
    for i in align_lines:
        if len(asequence) == len(bsequence):
            asequence.append(i)
        elif len(asequence) > len(bsequence):
            bsequence.append(i)
    #reformat asequence/bsequence into one giant string
    asequence = "".join(asequence)
    bsequence = "".join(bsequence)
    #make a pandas dataframe with one column for asequence, one column for b sequence, and one row for each character
    seq_align = pd.DataFrame({"asequence":list(asequence), "bsequence":list(bsequence)})
    #make a list starting at 1 with the first not dash in asequence, and increasing by 1 for each not dash in asequence, add - if asequence is a dash using pandas
    def counter_add(row, column_name):
        nonlocal counter
        if row[column_name] != "-":
            counter += 1
            return counter
        else:
            return counter
    counter = 0
    seq_align["asequence_index"] = seq_align.apply(counter_add, args=("asequence",), axis=1)
    counter = 0
    seq_align["bsequence_index"] = seq_align.apply(counter_add, args=("bsequence",), axis=1)
    return(seq_align)



#this function is getting replaced 
def mini_adjust_index(combined_list,index):
    for i in combined_list:
        if i[1] == index:
            return i[3]



def loop_flank2(res_info, pymoldict, pdbdict, refdict):
    #make a list of unique pdb_id and chain (pdb_chain)
    res_info=res_info.copy()
    #add 5 new columns to res_info: pdb_adj_res1, pdb_adj_res2, ref_pdb_chain, adj_res1, adj_res2
    res_info["pdb_adj_res1"] = ""
    res_info["pdb_adj_res2"] = ""
    res_info["ref_pdb_chain"] = ""
    res_info["adj_res1"] = ""
    res_info["adj_res2"] = ""
    res_info["flank1"] = ""
    res_info["flank2"] = ""
    unique_pdb_chain = res_info["pdb_chain"].unique()
    #align pymoldict and pdbdict with run_needle using pandas
    for i in unique_pdb_chain:
        #print progress with pdb_id and chain_id
        print(f"progress {unique_pdb_chain.tolist().index(i)+1}/{len(unique_pdb_chain)} {i}")
        #get pdb_id and chain_id
        pdb_id, chain_id = i.split("_")
        #get pymol sequence
        head1 = ">"+pdb_id + "_" + chain_id
        pymolsequence = pymoldict[head1]
        #get pdb sequence
        pdbsequence = pdbdict[head1]
        #run needle
        run_needle(pymolsequence,pdbsequence)
        #reformat alignment output
        align_df = alignment_format("align.txt")
        #delete align.txt
        os.remove("align.txt")
        #get pdb adjusted index for each res_index_res1 and res_index_res2 where res1_pdb_chain == i
        # Filter res_info to rows where res1_pdb_chain equals i
        filtered_res_info = res_info[res_info["pdb_chain"] == i].copy()
        #add new columns to filtered_res_info called pdb_adj_res1 and pdb_adj_res2
        filtered_res_info["pdb_adj_res1"] = ""
        filtered_res_info["pdb_adj_res2"] = ""
        #for each res_index_res1 and res_index_res2 in filtered_res_info, find the matching index in align_df
        #if align_df == "alignment failed", save "alignment failed" to pdb_adj_res1 and pdb_adj_res2
        if type(align_df) == str:
            res_info.loc[res_info["pdb_chain"] == i, "pdb_adj_res1"] = "alignment failed"
            res_info.loc[res_info["pdb_chain"] == i, "pdb_adj_res2"] = "alignment failed"
            res_info.loc[res_info["pdb_chain"] == i, "adj_res1"] = "alignment failed"
            res_info.loc[res_info["pdb_chain"] == i, "adj_res2"] = "alignment failed"
            #move on to the next index, row in res_info.iterros():
            continue
        for index, row in res_info.iterrows():
            if row["pdb_chain"] == i:
                #save res_index as variables
                res1 = row["res_index_res1"]
                res2 = row["res_index_res2"]
                #save adjusted index as variables
                res1_adj = align_df[align_df["asequence_index"] == res1]["bsequence_index"].values[0]
                res2_adj = align_df[align_df["asequence_index"] == res2]["bsequence_index"].values[0]
                #if res1_adj or res2_adj == 0, save "check alignment" to pdb_adj_res1 or pdb_adj_res2
                if res1_adj == 0 or res2_adj == 0:
                    res_info.at[index, "pdb_adj_res1"] = res1_adj
                    res_info.at[index, "pdb_adj_res2"] = res2_adj
                    res_info.at[index, "adj_res1"] = "check alignment"
                    res_info.at[index, "adj_res2"] = "check alignment"
                    #move on to the next index, row in res_info.iterros():
                    continue
                else:
                #add res1_adj as value in new column "pdb_adj_res1"
                    res_info.at[index, "pdb_adj_res1"] = res1_adj
                    res_info.at[index, "pdb_adj_res2"] = res2_adj
                # Check if matching indices are found
                flank1_start = max(0, int(res1_adj) - 6)
                flank1_end = min(len(pdbsequence), int(res1_adj) + 5)
                res_info.at[index, "flank1"] = pdbsequence[flank1_start:flank1_end]
            
                flank2_start = max(0, int(res2_adj) - 6)
                flank2_end = min(len(pdbsequence), int(res2_adj) + 5)
                res_info.at[index, "flank2"] = pdbsequence[flank2_start:flank2_end]
        #store the value of "mol_org" where pdb_chain == i as refdictstart
        refdictstart = filtered_res_info["mol_org"].values[0]
        #if align_df == "alignment failed", save "alignment failed" to pdb_adj_res1 and pdb_adj_res2
        #if align_df is a string, save "alignment failed" to adj_res1 and adj_res2
        if type(align_df) == str:
            res_info.loc[res_info["pdb_chain"] == i, "pdb_adj_res1"] = "alignment failed"
            res_info.loc[res_info["pdb_chain"] == i, "pdb_adj_res2"] = "alignment failed"
            res_info.loc[res_info["pdb_chain"] == i, "adj_res1"] = "alignment failed"
            res_info.loc[res_info["pdb_chain"] == i, "adj_res2"] = "alignment failed"
            #move on to the next index, row in res_info.iterros():
            continue
        #save new dictionary with only entrys where the key starts with refdictstart
        ref_key = ""
        refsequence = ""
        for key,value in refdict.items():
            if key.startswith(refdictstart):
                ref_key = key
                refsequence = value
                break
            else:
                ref_key = "no reference found"
                refsequence = "no reference found"
        #split ref_key by "|" and save the second element as ref_pdb_chain
        if ref_key == "no reference found":
            ref_pdb_chain = "no reference found"
        else:   
            ref_pdb_chain = ref_key.split("|")[1]
            res_info.loc[res_info["pdb_chain"] == i, "ref_pdb_chain"] = ref_pdb_chain
        #run needle
        run_needle(pdbsequence,refsequence)
        #reformat alignment output
        align_df = alignment_format("align.txt")
        #remove align.txt
        os.remove("align.txt")

        #get reference adjusted index for each res_index_res1 and res_index_res2 where res1_pdb_chain == i
        for index, row in res_info.iterrows():
            if row["adj_res1"] == "check alignment":
                continue
            if row["pdb_chain"] == i:
                #set ref_adj_res1 to the value in b_sequence index where asequence_index == pdb_adj_res1
                res_info.at[index, "adj_res1"] = align_df[align_df["asequence_index"] == row["pdb_adj_res1"]]["bsequence_index"].values[0]
                #set ref_adj_res2 to the value in b_sequence index where asequence_index == pdb_adj_res2
                res_info.at[index, "adj_res2"] = align_df[align_df["asequence_index"] == row["pdb_adj_res2"]]["bsequence_index"].values[0]
            
        #define res_info where pdb_chain == i for more info for trouble shooting
        #print(res_info[res_info["pdb_chain"] == i])
    return res_info

#make a new test_col equal to 1 if res_index_res1 is equal to adj_res1 and res_index_res2 is equal to adj_res2 using pandas
#res_info["test_col"] = res_info.apply(lambda row: 1 if row["res_index_res1"] == row["adj_res1"] and row["res_index_res2"] == row["adj_res2"] else 0, axis=1)        

        
final_resinfo = loop_flank2(df, pymol_dict, pdb_dict, reference_dict)


#write final_resinfo to csv
final_resinfo.to_csv(args.output + ".csv", index=False)








#~~~~~~~~~~~~~~~~ code for trouble shooting a broken function ~~~~~~~~~~~~~~~~~~~~~~~
def loop_flank_test(pdb_chain, res_info, pymoldict, pdbdict, refdict):
    #make a list of unique pdb_id and chain (pdb_chain)
    res_info=res_info.copy()
    #get_pdb_id and chain_id
    pdb_id, chain_id = pdb_chain.split("_")
    #get pymol sequence
    head1 = ">"+pdb_id + "_" + chain_id
    pymolsequence = pymoldict[head1]
    print(pymolsequence)
    #get pdb sequence
    pdbsequence = pdbdict[head1]
    print(pdbsequence)
    #run needle
    run_needle(pymolsequence,pdbsequence)
    #reformat alignment output
    align_df = alignment_format("align.txt")
    print(align_df)
    #get pdb adjusted index for each res_index_res1 and res_index_res2 where res1_pdb_chain == i
    # Filter res_info to rows where res1_pdb_chain equals i
    filtered_res_info = res_info[res_info["pdb_chain"] == pdb_chain].copy()
    #add new columns to filtered_res_info called pdb_adj_res1 and pdb_adj_res2
    filtered_res_info["pdb_adj_res1"] = ""
    filtered_res_info["pdb_adj_res2"] = ""
    #for each res_index_res1 and res_index_res2 in filtered_res_info, find the matching index in align_df
    for index, row in filtered_res_info.iterrows():
        if row["pdb_chain"] == pdb_chain:
            #save res_index as variables
            res1 = row["res_index_res1"]
            res2 = row["res_index_res2"]
            #save adjusted index as variables
            res1_adj = align_df[align_df["asequence_index"] == res1]["bsequence_index"].values[0]
            res2_adj = align_df[align_df["asequence_index"] == res2]["bsequence_index"].values[0]
            #if res1_adj or res2_adj == 0, save "check alignment" to pdb_adj_res1 or pdb_adj_res2
            if res1_adj == 0 or res2_adj == 0:
                filtered_res_info.at[index, "pdb_adj_res1"] = "check alignment"
                filtered_res_info.at[index, "pdb_adj_res2"] = "check alignment"
            else:
                #add res1_adj as value in new column "pdb_adj_res1"
                filtered_res_info.at[index, "pdb_adj_res1"] = res1_adj
                filtered_res_info.at[index, "pdb_adj_res2"] = res2_adj
    print("pdb_pymol alignment worked")
    #store the value of "mol_org" where pdb_chain == i as refdictstart
    refdictstart = filtered_res_info["mol_org"].values[0]
    #save new dictionary with only entrys where the key starts with refdictstart
    ref_key = ""
    refsequence = ""
    for key,value in refdict.items():
        if key.startswith(refdictstart):
            ref_key = key
            refsequence = value
            break
        else:
            ref_key = "no reference found"
            refsequence = "no reference found"
    #split ref_key by "|" and save the second element as ref_pdb_chain
    ref_pdb_chain = ref_key.split("|")[1]
    filtered_res_info.loc[filtered_res_info["pdb_chain"] == pdb_chain, "ref_pdb_chain"] = ref_pdb_chain
    #run needle
    run_needle(pdbsequence,refsequence)
    print("pdb_ref alignment worked")
    #reformat alignment output
    align_df = alignment_format("align.txt")
    print(align_df)
    print(filtered_res_info)
    #get reference adjusted index for each res_index_res1 and res_index_res2 where res1_pdb_chain == i
    for index, row in filtered_res_info.iterrows():
        if row["pdb_chain"] == pdb_chain:
            #if pdb_adj_res1 or pdb_adj_res2 == "check alignment", save "check alignment" to adj_res1 or adj_res2
            if row["pdb_adj_res1"] == "check alignment" or row["pdb_adj_res2"] == "check alignment":
                filtered_res_info.at[index, "adj_res1"] = "check alignment"
                filtered_res_info.at[index, "adj_res2"] = "check alignment"
            else:
                #set ref_adj_res1 to the value in b_sequence index where asequence_index == pdb_adj_res1
                filtered_res_info.at[index, "adj_res1"] = align_df[align_df["asequence_index"] == row["pdb_adj_res1"]]["bsequence_index"].values[0]
                #set ref_adj_res2 to the value in b_sequence index where asequence_index == pdb_adj_res2
                filtered_res_info.at[index, "adj_res2"] = align_df[align_df["asequence_index"] == row["pdb_adj_res2"]]["bsequence_index"].values[0]
        #define res_info where pdb_chain == i for more info for trouble shooting
        #print(res_info[res_info["pdb_chain"] == i])
    print(filtered_res_info)
  
#final_resinfo = loop_flank_test("5AFI_y",df, pymol_dict, pdb_dict, reference_dict)
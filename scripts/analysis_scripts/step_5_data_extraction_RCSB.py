#importing required libraries
import os
import pandas as pd
import numpy as np
import requests
from bs4 import BeautifulSoup
from optparse import OptionParser

#fixing '5E54' becoming '5.00E+54' issue
def fix_PDB_ID(df):
    for i, j in enumerate(df['PDB_ID']):
        print (j)
        if len(str(j))>4:
            if '-' in j:
                X= ''.join(j.split('-')).upper()
                df.loc[i, 'PDB_ID']= X
            else:
                if '.' in j:
                    x= str(j)
                    x1= x.split('.')[0]
                    x2= x.split('+')[1]
                    X= x1+'E'+x2
                    print (X)
                    df.loc[i, 'PDB_ID']= X
                else:
                    x= str(j)
                    x1= x.split('+')[0]
                    x2= x.split('+')[1]
                    X= x1+ x2
                    print (X)
                    df.loc[i, 'PDB_ID']= X
                    
    return (df)


#This function will extract detail information of each entity of corresponding PDB_ID from RCSB protein data bank
#last updated (February 2, 2024)
def get_info(p):
    en_dic={} #this dictionary will be populated for each entity id with chain list, molecule, organism, sequence, experiment, resolution information
    filename= p+'.cif'
    url= 'https://files.rcsb.org/view/%s' %filename #this url will extract experiment and resolution information 
    url_org = 'https://www.rcsb.org/structure/%s' %p #this url will extract experiment and resolution information
    url_seq= 'https://www.rcsb.org/fasta/entry/%s/display' %p #this url will extract information comes with fasta
    
    u = requests.get(url).content
    u_org = requests.get(url_org).content
    u_seq= requests.get(url_seq).content #fasta url for getting sequence information
    
    soup= BeautifulSoup(u_org, 'lxml')
    soup_seq= BeautifulSoup(u_seq, 'lxml')
    
    all_tables= pd.read_html(url_org) #in this list of tables, each tables contain detail information of each entity
    #however, we are using these tables only to extract chain ID and segment ID
    
    
    org= soup.find_all('li', id= 'header_organism') #html contents for source organism information
    se= soup_seq.find_all('body') #html contents for fasta

    #extracting experimental method information
    u1= str(u).split('#')
    for u2 in u1:  
        if '_exptl.method' in u2:
            exp_method=  u2.split("'")[1].replace('\\', '')
    
    #print (exp_method)
    #extracting resolution information 
    for u2 in u1:
        u3= u2.replace('_reconstruction.resolution_method', '')
        u4= u3.replace("'_refine.ls_d_res_high\\'", '')
        if '_em_3d_reconstruction.resolution' in u4:
            
            resolution1= (u4.replace('\\n', '').split('_em_3d_reconstruction.resolution')[1].split('_em_3d')[0]).replace(' ', '')
            if len(resolution1)==0:
                pass
            else:
                resolution= round(float(resolution1), 1)

        if '_refine.ls_d_res_high' in u4:
            print (u4)
            resolution1= (u4.replace('\\n', '').split('_refine.ls_d_res_high')[1].split('_refine.')[0]).replace(' ', '')
            if len(resolution1)==0:
                pass
            elif resolution1 =='.':
                pass
            else: 
                resolution= round(float(resolution1), 1)
            
         


    #getting detail information from fasta url
    for i in se:
        en_det=(str(i)[13: len(str(i))-11]).split('&gt;') #list of detail for each entity

    for i, j in enumerate(en_det):
        entity_id= (j.split('|')[0]).split('_')[1] #entity ID is a number, one for each entity.
        en_dic[entity_id]=[]
        #print (entity_id)
        
        for t1, t2 in enumerate(all_tables):
            if t2.shape== (1, 4): #protein entity table shape = (5, 6), nucleic acid entity table shape= (3, 6)
                pass
            else:
                #print (t2)
                if (type(t2.columns[0]))== tuple:
                    if t2.columns[0][0].split(':')[1][1:] == str(entity_id):
                        total_chain= (t2[t2.columns[1]][0])
        
        #extracting chain_list and seg_list
        chain_list=[]#list of chain ID for each entity ID, one entity can contain multiple chains
        #all the chain in one entity is same (not sure about this fact)
        seg_list=[] #list of seg ID for each entity ID, one entity can contain multiple segment
    
        
        #print (total_chain)
        #some examplex of total chain
        #PB [auth C0], VC [auth D0]
        #A, EB [auth FB]
        #A, C, E, G, IA, C, E, G, I, K, M, O Less --->example PDB: 4RMO
        
        #seg ID and chain ID can be interchanged in some fasta files 
        #therefore, seg ID and chain ID will be always extracted from the website
        total_chain_frags= total_chain.split(',')
        #print (total_chain_frags)
        for f1, f11 in enumerate(total_chain_frags): 
            f2= f11.replace('Less', '')
            #print (f2)
            if ']' in f2:
                list_str_total_chain= list(f2.split(']'))
                for h1, h2 in enumerate(list_str_total_chain):
                    list_str_total_chain[h1]= h2.replace(',', '')
                    list_str_total_chain[h1]= h2.replace(' ', '')
                if '' in list_str_total_chain:
                    list_str_total_chain.remove('')
                else:
                    pass
                
                for h1, h2 in enumerate(list_str_total_chain):
                    hh2= h2.split('[auth')
                    chain_list.append(hh2[1])
                    seg_list.append(hh2[0])
            else: 
                chain_list.append(f2.replace(' ', ''))
                seg_list.append(f2.replace(' ', ''))

        molecule= (j.split('|')[2]) #type of molecule
    
        organism= str(j.split('|')[3]).split('\n')[0] #this expressed organism
        if organism == 'null':
            organism = 'Not available'
        
        if len(org)==0:
            org_2= 'Not available'
        else:
            #find if multiple source organism is present
            org_1= (list(str(org[0]).split('>')))
            #print (org_1)
            org_2=''
            for or1, or2 in enumerate(org_1):
                if or2.endswith('</a'):
                    #print (or2)
                    org_2 += or2[:len(or2)-3]+ ', '
            org_2= org_2[: len(org_2)-2] #this is source organism
            #print (org_2)
    
        sequence= (j.split('|')[3]).split('\n')[1] #this the RNA sequence in expre
        #print (entity_id)
        #print (chain_list)
        #print (seg_list)
        
        en_dic[entity_id].append(chain_list) #chain list
        en_dic[entity_id].append(seg_list) #segment list
        en_dic[entity_id].append(molecule) #RNA type or protein type
        en_dic[entity_id].append(organism) #expressed organism
        en_dic[entity_id].append(org_2) #source organism
        en_dic[entity_id].append(sequence) #sequence
        en_dic[entity_id].append(exp_method) #experimental method
        en_dic[entity_id].append(resolution) #resolution
        #print (sequence)
    #print (en_dic)
    #en_dic[entity_id][0] --> chain list
    #en_dic[entity_id][1] --> segment list
    #en_dic[entity_id][2] --> molecule (eg. RNA type)
    #en_dic[entity_id][3] --> expressed organism
    #en_dic[entity_id][4] --> source organism
    #en_dic[entity_id][5] --> sequence
    #en_dic[entity_id][6] --> experimental method
    #en_dic[entity_id][7] --> resolution 
    
    return (en_dic)       

#this function will remove two types of problematic examples
#examples with different chain ID
#examples with icode in the residue index
def clean_two_prob(df):
    #converting residue index to string
    df['res_index_res1']= df["res_index_res1"].map(str)
    df['res_index_res2']= df["res_index_res2"].map(str)
    
    #Filter 1: remove examples where 'G' and 'U' are from different chain
    df1= df[df['chain_ID_res1'] == df['chain_ID_res2']]
    df1.index = np.arange(0, len(df1))
    
    #Filter 2: remove examples with icode
    for i, j in enumerate(df1['res_index_res1']):
        if '^' in j:
            df1.loc[i, 'res_index_res1']= '&&&'
        if '^' in df1['res_index_res2'][i]:
            df1.loc[i, 'res_index_res2']= '&&&'
    df2= df1[(df1['res_index_res1'] != '&&&') & (df1['res_index_res2'] != '&&&')]
    df2.index= np.arange(0, len(df2))
    
    return df2

# This function will create a dictionary where keys will be PC_ID (PDB ID_chain ID) and 
#values are list of ---> [0]chain ID list
#----------------------> [1]segmend ID list
#----------------------> [2]molecule (eg. RNA type)
#----------------------> [3]expressed organism
#----------------------> [4]source organism
#----------------------> [5]sequence
#----------------------> [6]experimental method
#----------------------> [7]resolution

def dict_pc_id(df):
    unq_pc_IDs= list(set(df['data_ext'].tolist()))

    print (unq_pc_IDs)
    D_pc_info={}
    for i, j in enumerate(unq_pc_IDs):
        print ('======================> working with this unique ID' + j)
        p= j.split('_')[0] #PDB ID
        c= j.split('_')[1] #Chain ID
        D_PDB_all= get_info(p)
        for k in D_PDB_all: #k is the entity ID
            if c in D_PDB_all[k][0]:
                print (j)
                D_pc_info[j]= D_PDB_all[k]

    return D_pc_info

#Step 1: preparing the dataframe
#read the data
optparser = OptionParser()
(options, args) = optparser.parse_args()
csvfile = args[0]

R= args[1]
D1= pd.read_csv(csvfile)

#we are checking now if there is any PDB_ID issue like '5E54' becoming '5.00E+54'
D2= fix_PDB_ID(D1)

#remove examples where residues are from different chain
#remove examples where residue index contain icode
D3= clean_two_prob(D2)

#droping chain_ID_res2 column
#because from now on all examples considered will contain residues from same chain
D4= D3.drop(['chain_ID_res2'], axis=1)
D4.rename(columns = {'chain_ID_res1':'chain_ID'}, inplace = True)

#inserting empty columns which will be populated in next steps
D4.insert (1, "Experimental_Method", '')
D4.insert (2, "Resolution_(Å)", '')
D4.insert(4, 'seg_ID', '')
D4.insert (5, 'Molecule', '')
D4.insert (6, 'Source_Organism', '')
D4.insert (7, 'Expressed_Organism', '')
D4.insert (8, 'Chain_length_reference', '')

D4.insert (9, 'data_ext', D4['PDB_ID']+'_'+D4['chain_ID']) #location of this column does not matter, because at the end we are going to delete this column


#Step 2: create a dictionary, where keys will be PC_ID (PDB ID_chain ID) from 'data_ext' column
# and value will be detail information for the corresponding PDB ID.
#these detail information will be extracted by get_info function
#in this entire script/notebook this is the most data intense step, therefore, will take the longest 
D_pc_info= dict_pc_id(D4)

#Step 3: lets populate those empty columns we have created on step 1
for i, j in enumerate(D4['data_ext']):
    print ('==========================================')
    print (j)
    chain_ID= j.split('_')[1] #this is chain ID
    print (chain_ID)
    print (D_pc_info)
    list_chain_ID1= (D_pc_info[j][0])
    list_chain_ID=[]
    for ch1, ch2 in enumerate(list_chain_ID1):
        ch3= ch2.replace(' ', '')
        list_chain_ID.append(ch3)
    print ('WORKING HERE!!!')
    print (list_chain_ID)
    chain_ID_index= list_chain_ID.index(chain_ID) #this is chain ID index in the D_pc_info, to find the segment ID the same index can be used
        
    D4.loc[i, 'Experimental_Method']= D_pc_info[j][6]
    D4.loc[i, 'Resolution_(Å)']= D_pc_info[j][7]  
    D4.loc[i, 'seg_ID']= D_pc_info[j][1][chain_ID_index]     
    D4.loc[i, 'Molecule']= D_pc_info[j][2]
    D4.loc[i, 'Source_Organism']= D_pc_info[j][4]
    D4.loc[i, 'Expressed_Organism']= D_pc_info[j][3]
    D4.loc[i, 'Chain_length_reference']= len(D_pc_info[j][5])


#Step 4: any additional data cleaning? 
#we will delete the 'data_ext' column now
D5= D4.drop(['data_ext'], axis=1)

#store the finalized data as csv
#store the finalized data as csv
if int(R)==1:
    #for standard wobble
    D5.to_csv('../../results/all_standard_wobble_RCSB_data.csv', index= False)
elif int(R)==2:
    #for shifted wobble
    D5.to_csv('../../results/all_shifted_wobble_RCSB_data.csv', index= False)
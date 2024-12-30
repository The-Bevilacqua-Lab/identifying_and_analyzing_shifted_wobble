#importing required libraries
import os
import pandas as pd
import numpy as np
import csv
import json
from os.path import isfile
from optparse import OptionParser

#required functions
#fixing '5E54' becoming '5.00E+54' issue
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

#directory of corresponding json files
#this dictionary will take the list of base pairs from the corresponding csv files
#and then will extract hydrogen bond information only for the base pairs of interest 
#for example, if we are interested on A-U base pairs, this function will only extract details
#of hydrogen bonds which are forming A-U base pairs

def pdb_hbond_dict(dir_j, test): 
    
    dir_w= os.getcwd()
    #this function will create a dictionary where:
    #keys are the PDB_ID
    #values are the dataframes which contain hydrogen bond details for the base pairs
    #these hydrogen bonds are the hydrogen bonds which help to form the base pairs of interest 
    #not all the hydrogen bonds
    
    #dir_j= directory of all json files
    #test= datafrmae of corresponding base pairs (G-C, A-U, G-U)
    
    unq_pdb= list(set(list(test['PDB_ID']))) #this is the unique PDB_IDs from the base pair dataframe
    pdb_bp={} #key= PDB_ID, value= base pair dataframe for the corresponding PDB_ID
    for k, l in enumerate(unq_pdb):
        pdb_bp[l]= test[test['PDB_ID']== l]
    print (pdb_bp)
    #reindexing dataframes inside pdb_bp
    for k1 in pdb_bp:
        pdb_bp[k1].index= np.arange(0, len(pdb_bp[k1]))
    #d3.index = np.arange(0, len(d3))
    
    os.chdir(dir_j)
    pdb_hbonds= {} #key= PDB_ID, value= hbond datafrme for hbonds forming base pairs

    for i, j in enumerate((unq_pdb)):
        print (i) #index
        print (j) #PDB_ID
        
        cif_n= j+'-dssr.json' #name of the json file for corresponding PDB_ID
        
        with open(cif_n, 'r') as f:
            data= json.loads(f.read())
            if 'pairs' in data:
                #the following dataframe contain all hbond information for the corresponding PDB_ID
                hb_list = pd.json_normalize(data, record_path =['hbonds'])
        
        #dictionary of bp names from different rules (name, Saenger, LW, DSSR)
        bp_name= {}
        for x1, x2 in enumerate(pdb_bp[j]['bp_res']):
            bp_name[x2]= []
            #name
            bp_name[x2].append(pdb_bp[j]['Name'][x1])
            #Saenger
            bp_name[x2].append(pdb_bp[j]['Saenger'][x1])
            #LW
            bp_name[x2].append(pdb_bp[j]['LW'][x1])
            #DSSR
            bp_name[x2].append(pdb_bp[j]['DSSR'][x1])
         
        
        
        #res1-->atom ID
        hb_list['atom_ID_res1']= hb_list['atom1_id'].str.split('@', expand = True)[0]
        #the 'extra_res1' column is just helping to create the 'chain ID', 'residue index', and 'residue ID' columns
        #this column will be deleted once all the column creation is done
        hb_list['extra_res1']= hb_list['atom1_id'].str.split('@', expand = True)[1]
        #res1-->chain ID
        hb_list['chain_ID_res1']= hb_list['extra_res1'].str.split('.', expand = True)[0]
        #res1-->residue index
        hb_list['res_index_res1']= hb_list['extra_res1'].str.split('.', expand = True)[1].str[1:]
        #res1-->residue ID
        hb_list['res_ID_res1']= hb_list['extra_res1'].str.split('.', expand = True)[1].str[0]
        
        #res2-->atom ID
        hb_list['atom_ID_res2']= hb_list['atom2_id'].str.split('@', expand = True)[0]
        #the 'extra_res2' column is just helping to create the 'chain ID', 'residue index', and 'residue ID' columns
        #this column will be deleted once all the column creation is done
        hb_list['extra_res2']= hb_list['atom2_id'].str.split('@', expand = True)[1]
        #res2-->chain ID
        hb_list['chain_ID_res2']= hb_list['extra_res2'].str.split('.', expand = True)[0]
        #res2-->residue index
        hb_list['res_index_res2']= hb_list['extra_res2'].str.split('.', expand = True)[1].str[1:]
        #res2-->residue ID
        hb_list['res_ID_res2']= hb_list['extra_res2'].str.split('.', expand = True)[1].str[0]
        
        
        #below we are adding one column 'bp_res', this column will help to exclude hbonds not part of the base pair of interest        
        #notation for this column: S.U625-S.G645
        hb_list['bp_res']= hb_list['atom1_id'].str.split('@',expand = True)[1]+'-' \
        +hb_list['atom2_id'].str.split('@',expand = True)[1]
        
        #Below we are adding one column 'bp_atom', this column will help to sort out different registers
        #notation for this column: G.N2-U.O4
        #we only need this column for G-U or U-G base pairs 
        #because we are not finding different registers from A-U, U-A, G-C, and C-G
        if test['Base_pair'][1] == 'A-C' or 'C-A':
            print ('#################################################')
            #print ((hb_list['atom1_id']))
            #print ((hb_list['atom1_id'].str.split('@', expand = True)[1].str.split('.', expand = True))[1][0])
            #print ((hb_list['atom1_id'].str.split('@', expand = True))[1])
            hb_list['bp_atom'] = hb_list['res_ID_res1']+'.'+hb_list['atom_ID_res1']+'-'+hb_list['res_ID_res2']+'.'+hb_list['atom_ID_res2']
        
        
        #below is the list which contains all the bps of the corresponding PDB_ID
        bps= list(pdb_bp[j]['bp_res']) 
        ##print (bps)
        
        #below we are only filtering out the bp_list by including only the hbonds involving the base pairs of interest
        hb_list_n = hb_list[hb_list['bp_res'].isin(bps)]
        hb_list_n.index = np.arange(0, len(hb_list_n))
        
        #now we are going to include the base pair name columns
        hb_list_n['name']=''
        hb_list_n['Saenger']= ''
        hb_list_n['LW']= ''
        hb_list_n['DSSR']= ''
        
        for x3, x4 in enumerate(hb_list_n['bp_res']):
            #name
            hb_list_n.loc[x3, 'name']= bp_name[x4][0]
            #Saenger
            hb_list_n.loc[x3, 'Saenger']= bp_name[x4][1]
            #LW
            hb_list_n.loc[x3, 'LW']= bp_name[x4][2]
            #DSSR
            hb_list_n.loc[x3, 'DSSR']= bp_name[x4][3]
        
        #adding a new column which will help to filter out registers
        hb_list_n.insert(0, 'PDB_ID', j)
        hb_list_n['bp_ID']= hb_list_n['PDB_ID'] +'_' +hb_list_n['bp_res']
        
        #deleting the 'extra_res1' and 'extra_res2' columns here
        hb_list_p= hb_list_n.drop(['extra_res1', 'extra_res2'], axis=1)
        
        
        
        #assigning the filtered hbond to the corresponding PDB_ID
        pdb_hbonds[j]= hb_list_p
    
    os.chdir(dir_w)
    return pdb_hbonds

#take the output from pdb_hbond dictionary and filter out different registers
#hbonds for different GU registers

def pdb_regs_dict(p):
    #hbonds for standard wobble
    #R1= ['G.O6-U.N3', 'G.N1-U.O2']
    #R2= ['U.N3-G.O6', 'U.O2-G.N1']

    #hbonds for shifted wobble
    R1= ['G.N1-U.O4', 'G.N2-U.N3']
    R2= ['U.O4-G.N1', 'U.N3-G.N2']

    
    r={} #dictionary contains hbonds for reg 2 
    for i in p:
        r[i]= p[i][(p[i]['bp_atom'].isin(R1))|(p[i]['bp_atom'].isin(R2))]
        r[i].index = np.arange(0, len(r[i]))

        
        #if no bp_atom is found in R21, R22, R31, and R32, it is possible to get empty dataframe
    
    R= [r]
    
    return R

#this function will take the pdb_regs dictionary
#this dictionary has PDB_ID as key and values will be register 2 oe register 3 examples within thecorresponding PDB_ID
#this function will concate the nonzero dataframes
def concat_regs(R):
    R_n=[]
    for i in (R):
        if len(R[i]) ==0:
            pass
        else:
            #R[i].insert(0, 'PDB_ID', i)
            RR= R[i].drop('index', axis= 1)
            R_n.append(RR)
    R_comb= pd.concat(R_n)
    
    return R_comb

def exclude_non_regs(d):
    #d is the dataframe contain all the base pairs which contain at least one of the hbonds 
    #present in rare wobble registers

    #for standard wobble
    R= ['G.O6-U.N3', 'G.N1-U.O2', 'U.N3-G.O6', 'U.O2-G.N1']
    #for shifted wobble
    R= ['G.N1-U.O4', 'G.N2-U.N3', 'U.O4-G.N1', 'U.N3-G.N2']
    #dataframe is grouped by one column (reg_2_ID) and then records from another column (bp_atom) for each group are 
    #stored as a list into a new column (list_hbond)
    d1= d.groupby('bp_ID')['bp_atom'].apply(list).reset_index(name= 'list_hbonds')
    
    D={} #in this dictionary key will be reg_2_ID and values will be the list of associated list of 'bp_atom'
    for i, j in enumerate(d1['bp_ID']):
        D[j]= d1['list_hbonds'][i]
        
        
    


    Rs=[] #reg_2_ID or reg_3_ID will be stored here if they contain the required hbonds for reg2 or reg3, respectively. 
    
    for i in D:
        if len(D[i])>1:
            print (i)
            check =  all(item in R for item in D[i])
            if check is True:
                print (D[i])
                Rs.append(i)
            else:
                print (D[i])
                pass

    d2= d[d['bp_ID'].isin(Rs)]
    d2.index = np.arange(0, len(d2))

    return d2

#read the csv file
optparser = OptionParser()
(options, args) = optparser.parse_args()
csvfile = args[0]
cif_dir = args[1]

D1= pd.read_csv(csvfile)

D2= fix_PDB_ID(D1)

D2['bp_res'] = D2['nt1']+'-'+D2['nt2']

D3= pdb_hbond_dict(cif_dir, D2)

D4= (pdb_regs_dict(D3))

D5= concat_regs(D4)

D6= exclude_non_regs(D5) #with hbonds

drp_cols= ['atom1_serNum', 'atom2_serNum', 'donAcc_type', 'distance', 'atom1_id', 'atom2_id', 'atom_pair', 'residue_pair', 'atom_ID_res1', 'atom_ID_res2', 'bp_atom']

D7= D6.drop_duplicates(subset=['bp_ID']) #drop duplicates 
D7.index= np.arange(0, len(D7)) #reindexing
D8= D7.drop(drp_cols, axis=1) #drop hbond columns 

#save output dataframes to data folder
#D6.to_csv('data/standard_GU_with_hbonds.csv', index= False)
#D8.to_csv('data/standard_GU_without_hbonds.csv', index=False)

#D6.to_csv('results/shifted_GU_with_hbonds.csv', index= False)
D8.to_csv('../../results/shifted_GU_without_hbonds.csv', index=False)
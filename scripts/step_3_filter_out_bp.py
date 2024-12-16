#importing required libraries
import os
import numpy as np
import pandas as pd
from os.path import isfile
from subprocess import call
from optparse import OptionParser
import json

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


optparser = OptionParser()
(options, args) = optparser.parse_args()
cif_dir = args[0]

D1= pd.read_csv('data/structures_within_3.2_resolution.csv')
print ('-----------------------------')
print (D1.shape)
D2= fix_PDB_ID(D1)


#STEP 1
#create a dictionary for PDB information
#for each PDB_ID as key, experimental method and resolution will be the value
pdb_info={}
for i, j in enumerate(D2['PDB_ID']):
    print (j)
    pdb_info[j] =[]
    pdb_info[j].append(D2['Experimental_Method'][i])
    pdb_info[j].append(D2['Resolution_(Å)'][i])
print (pdb_info)

#STEP 2
#this is the dataframe which will store all base pairs within all unique structure files
all_bps = pd.DataFrame(columns = ['PDB_ID', 'Experiment', 'Resolution', 'nt1', 'nt2', 'Base_pair','Name', 'Saenger', 'LW', 'DSSR'])
all_GUs = pd.DataFrame(columns = ['PDB_ID', 'Experiment', 'Resolution', 'nt1', 'nt2', 'Base_pair','Name', 'Saenger', 'LW', 'DSSR'])

#STEP 3
#we are changing the directory to cifs to open corresponding json files

home= os.getcwd() #this is our current working directory
os.chdir(cif_dir) #directory where the CIF and JSON files are stored

for ind, PID in enumerate(D2['PDB_ID']):
    print (ind)
    print (PID)
    cif_n= PID+'-dssr.json'
    #print (cif_n)
    with open(cif_n, 'r') as f:
        data= json.loads(f.read())
        #print (data)
        if 'pairs' in data:
            print ('it is working')
            bp_list = pd.json_normalize(data, record_path =['pairs'])

        if len(bp_list) ==0:
            pass
        else:
            for x1, x2 in enumerate(bp_list['bp']):
                #the next line is commented out now
                #running the next line can take very long time (weeks or months depending on the number of unique PDB_IDs)    
                all_bps.loc[len(all_bps)] = [j, pdb_info[j][0], pdb_info[j][1], bp_list['nt1'][x1], bp_list['nt2'][x1], bp_list['bp'][x1], bp_list['name'][x1], bp_list['Saenger'][x1], bp_list['LW'][x1], bp_list['DSSR'][x1]]
                
                if x2== 'G-U' or x2=='U-G': #this can be spacify as a variable before this for loop
                    all_GUs.loc[len(all_GUs)] = [j, pdb_info[j][0], pdb_info[j][1], bp_list['nt1'][x1], bp_list['nt2'][x1], bp_list['bp'][x1], bp_list['name'][x1], bp_list['Saenger'][x1], bp_list['LW'][x1], bp_list['DSSR'][x1]]
  
                    
os.chdir(home)

#print (all_bps)
#print (all_GUs)

all_GUs.to_csv('data/all_GU_base_pairs.csv', index= False)
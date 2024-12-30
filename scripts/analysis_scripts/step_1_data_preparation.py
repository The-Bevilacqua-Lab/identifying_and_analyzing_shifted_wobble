#importing required libraries
import os
import numpy as np
import pandas as pd
from optparse import OptionParser

#required functions
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
csv_dir = args[0]
R= float(args[1])

ds= []
#importing search data from RCSB-PDB
home= os.getcwd()
os.chdir(csv_dir)
for filename in os.listdir(csv_dir):
    if filename.endswith('.csv'):
        if filename.startswith('rcsb_pdb_custom_report_'):
            d= pd.read_csv(filename)
            ds.append(d)

os.chdir(home) 


# combining all search output df into one
D= pd.concat(ds)
D.index = np.arange(0, len(D))

D.columns = D.columns.str.replace(' ', '_') #replacing ' ' with '_' from all columns name

D= fix_PDB_ID(D)


#preparing a cleaner version of list of PDBs containing RNA molecules
D1 = D[['PDB_ID', 'Experimental_Method', 'Resolution_(Å)']].copy()

#because each rows in D are entities in the structures
#repeation of same structure (PDB_ID) is possible.
#now we are removing the duplicates, therefore, we will have detail information for each PDB_ID one time
D2= D1.drop_duplicates(subset=['PDB_ID'])
D2.index = np.arange(0, len(D2))

#filtering out empty cell at 'PDB_ID' columns
D3 = D2.dropna(subset=['PDB_ID'])

#sorting the dataframe by the PDB_ID
D4 = D3.sort_values(by='PDB_ID')
#print (D4.shape)


#filter out all base pair within a desired resolution cut-off
#apply resolution cut-off (3.2 A)
#removing any empty cells
D5= D4[~D4['Resolution_(Å)'].isna()]
D5.index = np.arange(0, len(D5))

#converting all entries in 'Resolution_(Å)' to float
D6= D5.copy()
D6['Resolution_(Å)'] = pd.to_numeric(D5['Resolution_(Å)'], errors='coerce')

#applying resolution cut-off
D7= D6[(D6['Resolution_(Å)']< float(R)) | (D6['Resolution_(Å)']== float(R)) ]
D7.index = np.arange(0, len(D7))

#D7 dataframe contain the PDB_IDs for the structures within user defined resolution cut-off
#this dataframe will be used in the next step to download and characterize the structures by DSSR
D7.to_csv('../../results/structures_within_'+ str(R)+ '_resolution.csv', index= False)
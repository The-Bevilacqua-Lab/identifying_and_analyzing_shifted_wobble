#importing required libraries
import os
import numpy as np
import pandas as pd
from os.path import isfile
from subprocess import call
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

#directory to store the cif files
#make sure this directory contain the executible DSSR
#read the csv file
optparser = OptionParser()
(options, args) = optparser.parse_args()
csvfile = args[0]
cif_dir = args[1]

D1= pd.read_csv(csvfile)
D2= fix_PDB_ID(D1)

#home directory where this script is stored
home= os.getcwd()

#read PDB_IDs from D7 and downloading in cif_dir
os.chdir(cif_dir)
for ind, PID in enumerate(D2['PDB_ID']):
    print (PID)
    
    pdbname= cif_dir+ '/'+ PID+ '.cif'
    if isfile(pdbname) == True:
        print ('FILE ALREADY EXISTS! NOT DOWNLOADING!')
    else: 
        print ('Downloading '+ PID)
        url = "http://www.rcsb.org/pdb/files/%s.cif" %PID
        call (["curl","-L","-O","-s",url])
    
os.chdir(home)

os.chdir(cif_dir)
stat_count= 0
for filename in os.listdir('.'):
    if filename.endswith('.cif'):
        stat_count +=1
        print (stat_count)
        pname= filename[:4]+'-dssr.json'
        x= cif_dir+'/'+filename[:4]+'-dssr.json'
        if isfile(x)== True:
            print ('ALREADY CHARACTERIZED')
        else:
            print ('NOT CHARACTERIZED YET')
            print (filename)
            os.system('./x3dna-dssr -i='+filename+ ' --json'+' -o='+ filename[:4]+'-dssr.json')
            for f in os.listdir('.'):
                if f.startswith('dssr-'):
                    os.remove(f)

os.chdir(home)
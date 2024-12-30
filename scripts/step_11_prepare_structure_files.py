#importing required libraries
import requests
import os
import pandas as pd
import numpy as np
import requests
from optparse import OptionParser

#functions
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

#function to extract structure file in cif format and convert as a pandas dataframe
def extract_cif(p):
    filename= p+'.cif'
    url= 'https://files.rcsb.org/view/%s' %filename 
    u = requests.get(url).content
    u1= str(u).split('#')
    if len(u1)==1:
        print ('STRUCTURE IS NOT AVAILABLE!')
        return '0'
    
    else:
        c_intrst= ['_atom_site.group_PDB', '_atom_site.id', '_atom_site.type_symbol',
           '_atom_site.label_atom_id', '_atom_site.label_alt_id',
           '_atom_site.label_comp_id', '_atom_site.label_asym_id',
           '_atom_site.label_entity_id', '_atom_site.label_seq_id',
           '_atom_site.pdbx_PDB_ins_code', '_atom_site.Cartn_x',
           '_atom_site.Cartn_y', '_atom_site.Cartn_z', '_atom_site.occupancy',
           '_atom_site.B_iso_or_equiv', '_atom_site.pdbx_formal_charge',
           '_atom_site.auth_seq_id', '_atom_site.auth_comp_id',
           '_atom_site.auth_asym_id', '_atom_site.auth_atom_id',
           '_atom_site.pdbx_PDB_model_num', '###']
        u2=[]
        for i, j in enumerate(u1):
            if c_intrst[0] in j:
                u2.append(j)

        if len(u2)==2:
            if len(u2[0])> len(u2[1]):
                u3= u2[0].replace('      ', '$')
            else:
                u3= u2[1].replace('      ', '$')
        else:
            u3= u2[0].replace('      ', '$')
        u4= u3.replace('     ', '$')
        u5= u4.replace('    ', '$')
        u6= u5.replace('   ', '$')
        u7= u6.replace('  ', '$')
        u8= u7.replace(' ', '$')
        u9= u8.replace('$$', '$')
        u10= u9.replace('\\\\', '').split('\\n')
        df = pd.DataFrame({'all_col':u10[23:]})

        df[c_intrst] = df['all_col'].str.split('$', expand=True)
        df= df.drop(['all_col', '###'], axis=1)

        return df
    
#function to clip residues of interest and store as cif file
def clip(p, c, s, l, n, h, m, F):
    #p= pdb_ID, entries can be case-insensitive
    #c= chain ID
    #s= seg ID, this is an optional value, if you do not know the seg ID of your residues of interest, just put ''
    #l= list of index for the residue of interest, items in this list should be strings, not integer
    #if you want all the residues for your chain of interest, just put '' for l
    #n= incase you have icode in the residue index for the residue of interest,
    #instead of providing all the residue indexes as a list in l, just give the index of the starting index
    #as a single string for l, and for n give the list of two numbers, first number will indicate how many residues 
    #before the residue of index you want, the second number will indicate how many residues after the residue of
    #inerest you want. If the residue of interest does not contain, icode, and/or you know all the residue indexes 
    #for the residues you want, just put '' for n
    #for instance, l= '47A', and n= [2, 7], means I want 2 residues before the residue at '47A' 
    # and 7 residue followed by residue at '47A'
    #h= 'N' or 'Y', if 'N', that means no hydrogen coordinates should be included considering hydrogen coordinates are provided
    #m= 'N' or 'Y', if 'N', that means no modified or non-RNA residues should be included 
    #if 'Y', that means no hydrogen coordinates should be removed if any hydrogen coordinates are provided
    #F= desired output filename (without the .cif extension)
    P= p.upper()
    
    #converting the cif file into dataframe
    p1= extract_cif(P)
    if len(p1)>1:
        #applying chain ID and seg ID filter 
        if len(s) ==0:
            p2= p1[(p1['_atom_site.auth_asym_id']==c)] 
            p2.index= np.arange(0, len(p2))
        else:
            p2= p1[(p1['_atom_site.auth_asym_id']==c)& (p1['_atom_site.label_asym_id']==s)] 
            p2.index= np.arange(0, len(p2))
    
        #converting index column values to integer
        #p2['_atom_site.auth_seq_id'] = pd.to_numeric(p2['_atom_site.auth_seq_id'])
        if isinstance(l, list):
            l= [str(item) for item in l]
            print (l)
            p3= p2[p2['_atom_site.auth_seq_id'].isin(l)]
            p3.index= np.arange(0, len(p3))
        elif len(l)>0:
            p2['index_icode']= p2['_atom_site.auth_seq_id']+ p2['_atom_site.pdbx_PDB_ins_code']
            #print (p2['index_icode'].to_list())
            ext_loc= int(p2[p2['index_icode'] == l]['_atom_site.label_seq_id'].to_list()[0])
            ext_locs= []
            for inds in range(ext_loc-n[0], ext_loc+n[1]):
                ext_locs.append(str(inds))
            p3= p2[p2['_atom_site.label_seq_id'].isin(ext_locs)]
            p3.index= np.arange(0, len(p3))
            p3 = p3.drop(columns=['index_icode'])
        elif len(l)==0:
            p3= p2.copy()
        
        #filtering out or not hydrogen atoms
        if h== 'N':
            p3_h= p3[p3['_atom_site.type_symbol']!= 'H']
            p3_h.index= np.arange(0, len(p3_h))
        else:
            p3_h= p3.copy()

        if m== 'N':
            p3_m = p3_h[p3_h['_atom_site.auth_comp_id'].str.len() < 2]
            p3_m.index= np.arange(0, len(p3_m))
        else:
            p3_m= p3_h.copy()
    
        #cleaning name for the sugar atoms
        for i, j in enumerate(p3_m['_atom_site.label_atom_id']):
            if "\\'" in j:
                #print (j)
                #print (j)
                x= (j.replace("\\'", "'"))
                y= (x.replace('"', ''))
                p3_m.loc[i, '_atom_site.label_atom_id'] = y
                p3_m.loc[i, '_atom_site.auth_atom_id'] = y
        
        #temp step
        #for now give everyting chain ID = 'A'
        p3_m['_atom_site.auth_asym_id']= 'A'
        #remove this section
    
        #Creating dataframe for the column name
        c_intrst_1= ['_atom_site.group_PDB', '_atom_site.id', '_atom_site.type_symbol',
               '_atom_site.label_atom_id', '_atom_site.label_alt_id',
               '_atom_site.label_comp_id', '_atom_site.label_asym_id',
               '_atom_site.label_entity_id', '_atom_site.label_seq_id',
               '_atom_site.pdbx_PDB_ins_code', '_atom_site.Cartn_x',
               '_atom_site.Cartn_y', '_atom_site.Cartn_z', '_atom_site.occupancy',
               '_atom_site.B_iso_or_equiv', '_atom_site.pdbx_formal_charge',
               '_atom_site.auth_seq_id', '_atom_site.auth_comp_id',
               '_atom_site.auth_asym_id', '_atom_site.auth_atom_id',
               '_atom_site.pdbx_PDB_model_num']
        c_intrst_2= ['data_', '#','loop_', '_atom_site.group_PDB', '_atom_site.id', '_atom_site.type_symbol',
               '_atom_site.label_atom_id', '_atom_site.label_alt_id',
               '_atom_site.label_comp_id', '_atom_site.label_asym_id',
               '_atom_site.label_entity_id', '_atom_site.label_seq_id',
               '_atom_site.pdbx_PDB_ins_code', '_atom_site.Cartn_x',
               '_atom_site.Cartn_y', '_atom_site.Cartn_z', '_atom_site.occupancy',
               '_atom_site.B_iso_or_equiv', '_atom_site.pdbx_formal_charge',
               '_atom_site.auth_seq_id', '_atom_site.auth_comp_id',
               '_atom_site.auth_asym_id', '_atom_site.auth_atom_id',
               '_atom_site.pdbx_PDB_model_num']
    
        c_intrst_2[0]= 'data_'+P
    
        p4= pd.DataFrame(columns= c_intrst_1)
        p4['_atom_site.group_PDB']= c_intrst_2 #this p3 dataframe above p2 will be joined in p4. p3 will help pymol to visualize the file
    
        p5 = pd.concat([p4, p3_m], axis=0)
    
        p5.to_csv(F+'.cif', header= False, sep=' ', index= False) #unmute this
        return p3_m
        #print ('WE COMPLETED THIS FARRRRRRRR')
        #p5.to_csv(P+'_clipped.cif', header= False, sep=' ', index= False)
    else:
        print ('STRUCTURE IS NOT AVAILABLE!')

def remove_mod_nts(tcif, p):
    #tcif is the structure file as a pandas dataframe format
    #p is the pdb ID
    
    #removing modified nucleotides
    cif = tcif[tcif['_atom_site.auth_comp_id'].str.len() <= 2]
    cif.index= np.arange(0, len(cif))
    
    #cleaning name for the sugar atoms
    for i, j in enumerate(cif['_atom_site.label_atom_id']):
        if "\\'" in j:
            #print (j)
            #print (j)
            x= (j.replace("\\'", "'"))
            y= (x.replace('"', ''))
            cif.loc[i, '_atom_site.label_atom_id'] = y
            cif.loc[i, '_atom_site.auth_atom_id'] = y
            
    #Creating dataframe for the column name
    c_intrst_1= ['_atom_site.group_PDB', '_atom_site.id', '_atom_site.type_symbol',
               '_atom_site.label_atom_id', '_atom_site.label_alt_id',
               '_atom_site.label_comp_id', '_atom_site.label_asym_id',
               '_atom_site.label_entity_id', '_atom_site.label_seq_id',
               '_atom_site.pdbx_PDB_ins_code', '_atom_site.Cartn_x',
               '_atom_site.Cartn_y', '_atom_site.Cartn_z', '_atom_site.occupancy',
               '_atom_site.B_iso_or_equiv', '_atom_site.pdbx_formal_charge',
               '_atom_site.auth_seq_id', '_atom_site.auth_comp_id',
               '_atom_site.auth_asym_id', '_atom_site.auth_atom_id',
               '_atom_site.pdbx_PDB_model_num']
    c_intrst_2= ['data_', '#','loop_', '_atom_site.group_PDB', '_atom_site.id', '_atom_site.type_symbol',
               '_atom_site.label_atom_id', '_atom_site.label_alt_id',
               '_atom_site.label_comp_id', '_atom_site.label_asym_id',
               '_atom_site.label_entity_id', '_atom_site.label_seq_id',
               '_atom_site.pdbx_PDB_ins_code', '_atom_site.Cartn_x',
               '_atom_site.Cartn_y', '_atom_site.Cartn_z', '_atom_site.occupancy',
               '_atom_site.B_iso_or_equiv', '_atom_site.pdbx_formal_charge',
               '_atom_site.auth_seq_id', '_atom_site.auth_comp_id',
               '_atom_site.auth_asym_id', '_atom_site.auth_atom_id',
               '_atom_site.pdbx_PDB_model_num']
    
    c_intrst_2[0]= 'data_'+ p
    
    cif1= pd.DataFrame(columns= c_intrst_1)
    cif1['_atom_site.group_PDB']= c_intrst_2
    
    cif2 = pd.concat([cif1, cif], axis=0)
    cif2.to_csv(p +'_mod_removed.cif', header= False, sep=' ', index= False)


#read the data
optparser = OptionParser()
(options, args) = optparser.parse_args()
csvfile = args[0]
D1= pd.read_csv(csvfile)
hom_dir= os.getcwd()
clipped_cif_dir= args[1] #this is the directory where you want to store the clipped and modified residue removed structure files 

D2= fix_PDB_ID(D1)

sw_df1= D2[D2['representative']==1]
sw_df1.index= np.arange(0, len(sw_df1))

unq_PDB_IDs= list(set(list(sw_df1['PDB_ID'].to_list()))) #list of unique pdb_ids
os.chdir(clipped_cif_dir)
for ind, pid in enumerate(unq_PDB_IDs):
    tcif= extract_cif(pid)
    tcif1= remove_mod_nts(tcif, pid)
os.chdir(hom_dir)
#importing required libraries
import os
import pandas as pd
import requests
import numpy as np
from Bio.PDB import MMCIFParser, Superimposer, PDBIO
from Bio.PDB.QCPSuperimposer import QCPSuperimposer
from optparse import OptionParser
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from collections import Counter

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
        
        #renumber the residue index starting from 1
        
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

def generate_str(df, tar_dir, fln):
    #df= dataframe for shifted wobble
    #tar_dir= directory where we want to store all clipped structures
    #fln= how many residues before and after the wobble residues should be included 

    home= os.getcwd() #home directory where this script is stored
    
    os.chdir(tar_dir)

    nr_df= df[df['representative']==1]
    nr_df.index= np.arange(0, len(nr_df))

    for k, l in enumerate(nr_df['bp_ID']):
        res1_index= int(nr_df['res_index_res1'][k])
        res2_index= int(nr_df['res_index_res2'][k])
        p= nr_df['PDB_ID'][k]
        c= nr_df['chain_ID'][k]
        s= nr_df['seg_ID'][k]
        resn_list_motif= [res1_index-fln, res1_index, res1_index+fln, res2_index-fln, res2_index, res2_index+fln]
        resn_list_motif.sort()
        filename= l
        cls_str= clip(p, c, s, resn_list_motif, '', 'N', 'Y', filename)

    os.chdir(home)


##############THE FUNCTIONS ABOVE ARE FOR GENERATING THE STRUCTURES, WHICH WILL ALIGNED AGAINST EACH OTHER WITH THE FUNCTIONS BELOW

def who_is_missing(structure, r1, r2):
    #r1 is the first residue of wobble
    #r2 is the second residue of wobble
    dict_good_length={'A':22, 'G':23, 'OMG':23, 'C':20, 'U':20, 'OMU':20}
    residues_s= []
    residues_e= {r1-1: 1, r1: 2, r1+1: 3, r2-1: 4, r2: 5, r2+1: 6}
    res_index_s=[]
    for model in structure:
        for chain in model:
            for residue in chain:
                atoms= list(residue.get_atoms())
                atoms_noH=[]
                #print (atoms)
                for atom in atoms:
                    if 'H' in atom.get_name():
                        pass
                    else:
                        atoms_noH.append(atom)
                if len(atoms_noH)== dict_good_length[residue.get_resname()]: #no atoms are missing, take this residue
                    residues_s.append(residue.get_id()[1])
                else: #some atom/atoms are missing, exclude this residue
                    pass

    #print (residues_s)
    for residue in residues_e:
        if residue in residues_s:
            res_index_s.append(residues_e[residue])
        else:
            pass
    #res_index_s contain the index for which nothing is missing
    return res_index_s

##residues1_not_missing= who_is_missing(structure1, 1086, 1099)
##residues2_not_missing= who_is_missing(structure2, 1083, 1096)


##common_index= list(set(residues1_not_missing).intersection(residues2_not_missing))
##print (common_index)
def req_atoms_alignment(structure, co_ind, res1, res2):
    atom_dicts= {'A': ["p", "C4'", "N9", "C6", "C2"], 'G': ["p", "C4'", "N9", "C6", "C2"], 'OMG': ["p", "C4'", "N9", "C6", "C2"],'C': ["p", "C4'", "N1", "C2", "C4"],'U': ["p", "C4'", "N1", "C2", "C4"], 'OMU': ["p", "C4'", "N1", "C2", "C4"]}
    residues_e= {res1-1: 1, res1: 2, res1+1: 3, res2-1: 4, res2: 5, res2+1: 6}
    count=0
    good_atoms=[]
    for model in structure:
        for chain in model:
            for residue in chain:
                if residues_e[residue.get_id()[1]] in co_ind:
                    atoms= list(residue.get_atoms())
                    #atoms_noH=[]
                    for atom in atoms:
                        #print ('----------------------------->.......look here')
                        #print (atom.get_name())
                        if atom.get_name() in atom_dicts[residue.get_resname()]:
                            good_atoms.append(atom)
                else:
                    pass
    
    return good_atoms


def align(str1, str2):
    parser = MMCIFParser(QUIET=True)
    superimposer = Superimposer()
    structure1 = parser.get_structure("RNA1", str1)
    structure2 = parser.get_structure("RNA2", str2)

    atom_dicts= {'A': ["p", "C4'", "N9", "C6", "C2"], 'G': ["p", "C4'", "N9", "C6", "C2"], 'C': ["p", "C4'", "N1", "C2", "C4"], 'U': ["p", "C4'", "N1", "C2", "C4"]}
    #reg3_8EV3_6.G105-6.U179.cif

    str1_res1= int(str1.split('_')[2].split('-')[0][1:].rstrip('.cif'))
    str1_res2= int(str1.split('_')[3][1:].rstrip('.cif'))

    str2_res1= int(str2.split('_')[2].split('-')[0][1:].rstrip('.cif'))
    str2_res2= int(str2.split('_')[3][1:].rstrip('.cif'))

    #find residues which will be considered to calculate RMSD
    str1_good_res= who_is_missing(structure1, str1_res1, str1_res2) #good residues in structure1
    str2_good_res= who_is_missing(structure2, str2_res1, str2_res2) #good residues in structure2

    common_goods= list(set(str1_good_res).intersection(str2_good_res)) #common

    structure1_atoms= req_atoms_alignment(structure1, common_goods, str1_res1, str1_res2)
    structure2_atoms= req_atoms_alignment(structure2, common_goods, str2_res1, str2_res2)

    superimposer.set_atoms(structure1_atoms, structure2_atoms)
    rmsd = superimposer.rms

    return rmsd


def align_all_against_all(dr):
    home= os.getcwd()
    os.chdir(dr)
    clmns=['#']
    filenames=[]
    RMSDs={}
    for filename in os.listdir('.'):
        if filename.endswith('.cif'):
            clmns.append(filename)
            filenames.append(filename)
    print (filenames)
        
    t1= pd.DataFrame(columns = clmns)
    t1['#']= filenames

    def unique_pairs(items):
        pairs = []
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                pairs.append(str(items[i])+'*'+str(items[j]))
        return pairs
    unq_pairs= unique_pairs(filenames)

    for i, j in enumerate(unq_pairs):
        #unq_pairs_RMSD[j]=[]
        k1= j.split('*')[0]
        k2= j.split('*')[1]

        t1.loc[t1['#']== k1, k2]= round(align(k1, k2), 2)
        t1.loc[t1['#']== k2, k1]= round(align(k2, k1), 2)
        t1.loc[t1['#']== k1, k1]= 0
        t1.loc[t1['#']== k2, k2]= 0
        if round(align(k1, k2), 2) != round(align(k2, k1), 2):
            print ('<><><><><><><><><><><><><><><><><><><><><><><>')
            print (j)
            print (align(k1, k2))
            print (align(k2, k1))

    
    os.chdir(home)
    return t1

#read the data
optparser = OptionParser()
(options, args) = optparser.parse_args()
csvfile = args[0]
D1= pd.read_csv(csvfile)

tar_dir= args[1]

fln= int(args[2])

#generating structures for alignment
generate_str(D1, tar_dir, fln)

#perform all-against-all alignment
D2= align_all_against_all(tar_dir)

D2.to_csv('results/all_against_all_RMSD_nr_shifted_wobbles.csv', index= False)

sw_names= [j.rstrip('.cif').replace('-', '_') for i,j in enumerate(list(D2.columns)[1:])]

D3= D2.drop(['#'], axis=1)
D3.index= sw_names
D3.columns= sw_names

# Convert DataFrame to squareform distance matrix
rmsd_matrix = squareform(D3.values)

# Perform hierarchical clustering
Z = linkage(rmsd_matrix, method='average')
Z[:, [0, 1]] = Z[:, [1, 0]] 

dend = dendrogram(Z, color_threshold= 1.23, labels=D3.index, above_threshold_color='lightgray', leaf_rotation=0, orientation='left', no_plot=True)


color_list_counter=  dict(Counter(dend['color_list']))
lcolor_list_counter=  dict(Counter(dend['leaves_color_list']))

small_clusters=[] #clusters with less than 4 members will be stored here

for i in lcolor_list_counter:
    if i != 'lightgray' and lcolor_list_counter[i]<4:
        small_clusters.append(i)
        
color_list_n= [] #updated color list for the branches
leaves_color_list_n= [] #updated color list for the branch leaves

for i, cls in enumerate(dend['color_list']):
    if cls in small_clusters:
        color_list_n.append('lightgray')
    else:
        color_list_n.append(cls)

for i, ls in enumerate(dend['leaves_color_list']):
    if ls in small_clusters:
      leaves_color_list_n.append('lightgray')
    else:
        leaves_color_list_n.append(ls)

dend["color_list"]= color_list_n #updated branch color list is assigned to the original dendrogram
dend["leaves_color_list"]= leaves_color_list_n #updated leaves color list is assigned to the original dendrogram

cl_dict={}

for cls1, cls2 in enumerate(dend['leaves_color_list']):
    if cls2== 'lightgray':
        pass
    else:
        if cls2 in cl_dict:
            cl_dict[cls2].append(dend['ivl'][cls1])
        else:
            cl_dict[cls2]= []
            cl_dict[cls2].append(dend['ivl'][cls1])      
print (cl_dict)

cl_dict= {f"{k[:1]}{i+1}": v for i, (k, v) in enumerate(cl_dict.items()) if k != 'lightgray'}
print ('-------------------------------------------------')
print (cl_dict)

def get_key(value):
    for key, values in cl_dict.items():
        if value in values:
            return key
D1['bp_ID']= D1['bp_ID'].str.replace('-', '_', regex=False)        
D1['cluster'] = D1['bp_ID'].apply(get_key)
D1['cluster'] = D1['cluster'].fillna('unclustered')

D1.to_csv('../results/nr_shifted_wobble_structural_clusters.csv', index= False)
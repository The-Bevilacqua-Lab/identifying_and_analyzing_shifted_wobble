import requests
import os
from subprocess import call
import pandas as pd
import csv
import json
import glob
import subprocess
import numpy as np
from Bio.PDB import MMCIFParser, Superimposer, PDBIO
from Bio.PDB.Structure import Structure
import Bio.PDB
from Bio.PDB.mmcifio import MMCIFIO
import shutil
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
        
        if len(u10[24].split('$')) == 22:
            df = pd.DataFrame({'all_col':u10[23:]})
            df[c_intrst] = df['all_col'].str.split('$', expand=True)
            df= df.drop(['all_col', '###'], axis=1)
            
        elif len(u10[24].split('$')) == 23:
            c_intrst[21]= '_atom_site.calc_flag'
            c_intrst.append('###')
            df = pd.DataFrame({'all_col':u10[24:]})
            df[c_intrst] = df['all_col'].str.split('$', expand=True)
            df= df.drop(['all_col', '###'], axis=1)

        return df

#function to clip residues of interest and store as cif file
def clip(p, c, s, l, n, h, F):
    #p= pdb_ID, entries can be case-insensitive
    #c= chain ID
    #s= seg ID, this is an optional value, if you do not know the seg ID of your residues of interest, just put ''
    #l= list of index for the residue of interest, items in this list should be strings, not integer
    #n= incase you have icode in the residue index for the residue of interest,
    #instead of providing all the residue indexes as a list in l, just give the index of the starting index
    #as a single string for l, and for n give the list of two numbers, first number will indicate how many residues 
    #before the residue of index you want, the second number will indicate how many residues after the residue of
    #inerest you want. If the residue of interest does not contain, icode, and/or you know all the residue indexes 
    #for the residues you want, just put '' for n
    #for instance, l= '47A', and n= [2, 7], means I want 2 residues before the residue at '47A' 
    # and 7 residue followed by residue at '47A'
    #h= 'N' or 'Y', if 'N', that means no hydrogen coordinates should be included considering hydrogen coordinates are provided
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
        else:
            p2['index_icode']= p2['_atom_site.auth_seq_id']+ p2['_atom_site.pdbx_PDB_ins_code']
            #print (p2['index_icode'].to_list())
            ext_loc= int(p2[p2['index_icode'] == l]['_atom_site.label_seq_id'].to_list()[0])
            ext_locs= []
            for inds in range(ext_loc-n[0], ext_loc+n[1]):
                ext_locs.append(str(inds))
            p3= p2[p2['_atom_site.label_seq_id'].isin(ext_locs)]
            p3.index= np.arange(0, len(p3))
            p3 = p3.drop(columns=['index_icode'])
        
        #filtering out or not hydrogen atoms
        if h=='N':
            p3_h= p3[p3['_atom_site.type_symbol']!= 'H']
            p3_h.index= np.arange(0, len(p3_h))
        else:
            p3_h= p3.copy()
    
        #cleaning name for the sugar atoms
        for i, j in enumerate(p3_h['_atom_site.label_atom_id']):
            if "\\'" in j:
                #print (j)
                #print (j)
                x= (j.replace("\\'", "'"))
                y= (x.replace('"', ''))
                p3_h.loc[i, '_atom_site.label_atom_id'] = y
                p3_h.loc[i, '_atom_site.auth_atom_id'] = y
    
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
        
        p_col= p3_h.columns
        if len(p_col)==22:
            c_intrst_1.append('_atom_site.calc_flag')
            c_intrst_2.append('_atom_site.calc_flag')
        else:
            pass
        
        c_intrst_2[0]= 'data_'+P
    
        p4= pd.DataFrame(columns= c_intrst_1)
        p4['_atom_site.group_PDB']= c_intrst_2 #this p3 dataframe above p2 will be joined in p4. p3 will help pymol to visualize the file
    
        p5 = pd.concat([p4, p3_h], axis=0)
    
        p5.to_csv(F+'.cif', header= False, sep=' ', index= False) #unmute this
        return p3_h
        print ('WE COMPLETED THIS FARRRRRRRR')
        #p5.to_csv(P+'_clipped.cif', header= False, sep=' ', index= False)
    else:
        print ('STRUCTURE IS NOT AVAILABLE!')

def load_structures_from_directory(directory, f_list):
    parser = MMCIFParser(QUIET=True)
    structures = []
    #file_paths = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.cif')]
    file_paths = [os.path.join(directory, f) for f in os.listdir(directory) if f in(f_list)]
    print (file_paths)

    for file_path in file_paths:
        structure = parser.get_structure(file_path, file_path)
        structures.append(structure)
    
    return structures


#function to check missing atoms of a clipped RNA motif
def find_miss_atoms(str_df):
    good_length= {'A': 22, 'C': 20, 'G': 23, 'U': 20, 'OMU':20, 'OMG': 24, '2MG': 24, '5MC': 21, '4AC': 23, 'PSU': 20, 'M7A': 23, 'A2M': 22, 'C2M': 20, 'G2M': 23, 'U2M': 20,}
    missing_atoms= []
    residues= list(set(str_df['_atom_site.auth_seq_id'].to_list()))
    print (residues)
    for i, j in enumerate(residues):
        res_atoms= str_df[str_df['_atom_site.auth_seq_id']== j]
        res_atoms.index= np.arange(0, len(res_atoms))
        res_ID= list(set(res_atoms['_atom_site.auth_comp_id'].to_list()))[0]
        if good_length[res_ID]== len(res_atoms):
            print ('NOTHING IS MISSING')
        else:
            missing_atoms.append(j)
    return missing_atoms   


def req_atoms_alignment(structure):
    atom_dicts= {'A': ["p", "C4'", "N9", "C6", "C2"], 'G': ["p", "C4'", "N9", "C6", "C2"], 'OMG': ["p", "C4'", "N9", "C6", "C2"], '2MG': ["p", "C4'", "N9", "C6", "C2"], 'C': ["p", "C4'", "N1", "C2", "C4"],'U': ["p", "C4'", "N1", "C2", "C4"], 'OMU': ["p", "C4'", "N1", "C2", "C4"], '5MC': ["p", "C4'", "N1", "C2", "C4"],'4AC': ["p", "C4'", "N1", "C2", "C4"], 'PSU': ["p", "C4'", "N1", "C2", "C4"], 'M7A': ["p", "C4'", "N9", "C6", "C2"], 'A2M': ["p", "C4'", "N9", "C6", "C2"], 'G2M': ["p", "C4'", "N9", "C6", "C2"], 'C2M': ["p", "C4'", "N1", "C2", "C4"],'U2M': ["p", "C4'", "N1", "C2", "C4"]}
    count=0
    good_atoms=[]
    for model in structure:
        for chain in model:
            for residue in chain:
                atoms= list(residue.get_atoms())
                #atoms_noH=[]
                for atom in atoms:
                    #print ('----------------------------->.......look here')
                    #print (atom.get_name())
                    if atom.get_name() in atom_dicts[residue.get_resname()]:
                        good_atoms.append(atom)
    
    print (good_atoms)
    return good_atoms


def align_structures(structures):
    print (structures)
    ref_structure = structures[0]
    #ref_atoms = [atom for atom in ref_structure.get_atoms()]
    ref_atoms = req_atoms_alignment(ref_structure)
    
    super_imposers = []
    
    for structure in structures[1:]:
        filename = os.path.basename(str(structure)).rstrip('.cif>')
        print(filename)
        #atoms = [atom for atom in structure.get_atoms()]
        atoms = req_atoms_alignment(structure)
        super_imposer = Superimposer()
        super_imposer.set_atoms(ref_atoms, atoms)
        #super_imposer.apply(structure.get_atoms().get_coords())
        
        super_imposer.apply(structure.get_atoms())

        # Save the aligned version of structures
        ##io = io=MMCIFIO()
        ##io.set_structure(structure) 
        ##io.save(filename+"_aligned.cif")
        
        print (super_imposer.rms)
        super_imposers.append(structure)
    
    return super_imposers


def compute_average_structure(structures, super_imposers):
    if not structures:
        raise ValueError("No structures provided.")
    
    ref_structure = structures[0]
    ref_atoms = req_atoms_alignment(ref_structure)
    num_atoms = len(ref_atoms)
    print ('look here')
    print (ref_atoms)
    
    coords_sum = np.zeros((num_atoms, 3))
    
    for i, structure in enumerate(structures):
        print (structure)
        if i > 0:
            #super_imposers[i-1].apply(structure.get_atoms().get_coords())
            #super_imposers[i-1].get_atoms().get_coords()
            atoms = [atom for atom in req_atoms_alignment(super_imposers[i-1])]
            coords = np.array([atom.get_coord() for atom in atoms])
            print (len(coords))
            coords_sum += coords
        else:
            atoms = [atom for atom in req_atoms_alignment(structure)]
            coords = np.array([atom.get_coord() for atom in atoms])
            print (len(coords))
            coords_sum += coords
    
    avg_coords = coords_sum / len(structures)
    #print ((avg_coords[23]))
    
    avg_structure = ref_structure.copy()

    for i, atom in enumerate(req_atoms_alignment(avg_structure)):
        print (i)
        print (atom)
        print ('----------------before')
        print (atom.get_coord())
        atom.set_coord(avg_coords[i])
        print ('----------------after')
        print (atom.get_coord())

    return avg_structure

def align_avg(avg, strs):
    avg_atoms= req_atoms_alignment(avg)
    rmsd_dict={}
    for structure in strs:
        filename = os.path.basename(str(structure)).rstrip('.cif>')
        print(filename)
        
        atoms = req_atoms_alignment(structure)
        super_imposer = Superimposer()
        super_imposer.set_atoms(avg_atoms, atoms)
        
        rmsd_dict[filename]= super_imposer.rms
        
    return rmsd_dict

#This function will remove the redundant data points
def rem_red(dfr):
    df= fix_PDB_ID(dfr)
    #---------------------------------------------
    for i, j in enumerate(df['adj_res1']):
        if pd.isnull(df.loc[i, 'adj_res1']):
            df.loc[i, 'adj_res1'] = df['res_index_res1'][i]
            df.loc[i, 'adj_res2'] = df['res_index_res2'][i]
        if df['adj_res1'][i]== df['adj_res2'][i]:
            df.loc[i, 'adj_res1'] = df['res_index_res1'][i]
            df.loc[i, 'adj_res2'] = df['res_index_res2'][i]

    df['adj_res1']= df['adj_res1'].map(str)
    df['adj_res2']= df['adj_res2'].map(str)

    df.insert (13, 'bp_notation', df['res_ID_res1']+ df['adj_res1'] +'-' +df['res_ID_res2']+ df['adj_res2'])
    #---------------------------------------------
    
    
    
    dfs=[]
    df_n= df[(df['b_resG']<=100) & (df['b_resU']<=100)]
    df_n.index= np.arange(0, len(df_n))
    ##print (df_n.shape)
    df1=df_n.groupby(['Molecule', 'Source_Organism_chain', 'bp_notation'])['Chain_length_structure'].describe().reset_index()
    ##print (df1.shape)
    
    red_ID=0
    work_stat=0
    for i, j in enumerate(df1['bp_notation']):
        work_stat +=1
        print (work_stat)
        print ('-----------------------------------------------------------------------------> working here '+ j)
        
        M= df1['Molecule'][i]
        S= df1['Source_Organism_chain'][i]
        B= df1['bp_notation'][i]
        df2= df_n[(df_n['Molecule']==M) & (df_n['Source_Organism_chain']==S) &(df_n['bp_notation']==B)]
        df2.index = np.arange(0, len(df2)) 
        red_ID +=1
        df2['red_group_ID']= red_ID
        df2['number_of_group_members']= len(df2)
        df2['motif_length']= '' #should be equal to 6 if no residue is missing
        df2['missing_res_atoms']=''
        df2['RMSD_to_avg']= ''
        df2['representative']= ''
        
        #0= nothing is missing
        #1= atoms missing
        #2= one entire residue is missing
        #3= both one or more entire residues and some atoms are missing
        
        #check missing residues for all examples
        #######################################################################
        for i, j in enumerate(df2['bp_ID']):
            p= df2['PDB_ID'][i]
            c= str(df2['chain_ID'][i])
            s= str(df2['seg_ID'][i])
            ind1= int(df2['res_index_res1'][i])
            ind2= int(df2['res_index_res2'][i])
            l= [ind1-1, ind1, ind1+1, ind2-1, ind2, ind2+1]
            n= ''
            h= 'N'
            F= j
            missing_residues= []
            #print (p)
            #print (c)
            #print (s)
            #print (l)
            #print (h)
            #print (F)
            print ('generating this structure: '+F)
            #print (p)
            str_file= clip(p, c, s, l, n, h, F)
            #print (len(str_file))

            df2.loc[(df2['bp_ID']==j), "motif_length"]= len(set(str_file['_atom_site.auth_seq_id'].to_list()))    
            if len(set(str_file['_atom_site.auth_seq_id'].to_list())) <6:
                print ('#####################################################')
                print (j)
                print ('#####################################################')
                #missing_residues[j]= len(set(str_file['_atom_site.auth_seq_id'].to_list()))
                missing_residues.append(j)
            else:
                pass
                #print (len(set(str_file['_atom_site.auth_seq_id'].to_list())))
            missing_atoms= find_miss_atoms(str_file)
            
            if len(missing_atoms)==0 and len(missing_residues)==0:
                df2.loc[(df2['bp_ID']==j), "missing_res_atoms"]= 0
                
            elif len(missing_atoms)>=1 and len(missing_residues)==0:
                df2.loc[(df2['bp_ID']==j), "missing_res_atoms"]= 1
            
            elif len(missing_atoms)==0 and len(missing_residues)==1:
                df2.loc[(df2['bp_ID']==j), "missing_res_atoms"]= 2
            
            elif len(missing_atoms)>=1 and len(missing_residues)==1:
                df2.loc[(df2['bp_ID']==j), "missing_res_atoms"]= 3
        
        #######################################################################
        if len(df2)== 1:
            df2['representative']= 1
            #df2.loc[((df2['bp_ID']==p1)&(df['bp_res']==x2)), "N2_N3_dih"]= N2_N3_dih
            dfs.append(df2)
            
        else:
            
            df_g= df2[df2["missing_res_atoms"]==0]
            df_g.index = np.arange(0, len(df_g)) 
            
            if len(df_g)==1:
                rep_bp_ID= df_g['bp_ID'][0]
                df2.loc[(df2['bp_ID']== rep_bp_ID), "representative"]= 1
                dfs.append(df2)
                
            elif len(df_g)==2:   
                b_list= df_g['b_comb'].to_list()
                if len(set(b_list))==1:
                    print (df_g['bp_ID'])
                    df3 = df_g.sample(n=1)
                    df3.index = np.arange(0, len(df3)) 
                else:
                    df3= df_g[df_g['b_comb']== min(b_list)]
                    df3.index = np.arange(0, len(df3)) 
            
                rep_bp_ID= df3['bp_ID'][0]
                df2.loc[(df2['bp_ID']== rep_bp_ID), "representative"]= 1
                dfs.append(df2)
                
            elif len(df_g)>2:

                str_dir= os.getcwd()
                str_list= [fname+ '.cif' for fname in df_g['bp_ID'].to_list()]
                
                structures1= load_structures_from_directory(str_dir, str_list) #loading required structures
                
                structures2= align_structures(structures1) #aligning structures with the reference structure
                
                avg_str= compute_average_structure(structures1, structures2) #calculating average structure
                
                rmsds= align_avg(avg_str, structures1) #calculating distance from the average structure
                
                rep_strs = min(rmsds, key=lambda k: rmsds[k])
                
                for bp_IDs in rmsds:
                    df2.loc[(df2['bp_ID']== bp_IDs), "RMSD_to_avg"]= rmsds[bp_IDs] #check this line again
                
                df2.loc[(df2['bp_ID']== rep_strs), "representative"]= 1
                dfs.append(df2)
                
                
            elif len(df_g)==0:
                #work with df2
                #none of the motifs are good
                motif_lengths= df2['motif_length'].to_list()
                df_o= df2[df2['motif_length']== max(motif_lengths)] #df okay, everyone is missing something, here we have the ones with minimum loss
                df_o.index= np.arange(0, len(df_o))

                if len(df_o)==1:
                    rep_bp_ID= df_o['bp_ID'][0]
                    df2.loc[(df2['bp_ID']== rep_bp_ID), "representative"]= 1
                    dfs.append(df2)
                elif len(df_o)==2:   
                    b_list= df_o['b_comb'].to_list()
                    if len(set(b_list))==1:
                        print (df_o['bp_ID'])
                        df3 = df_o.sample(n=1)
                        df3.index = np.arange(0, len(df3)) 
                    else:
                        df3= df_o[df_o['b_comb']== min(b_list)]
                        df3.index = np.arange(0, len(df3)) 
            
                    rep_bp_ID= df3['bp_ID'][0]
                    df2.loc[(df2['bp_ID']== rep_bp_ID), "representative"]= 1
                    dfs.append(df2)
                
                elif len(df_o)>2:
                    str_dir= os.getcwd()
                    str_list= [fname+ '.cif' for fname in df2['bp_ID'].to_list()]
                
                    structures1= load_structures_from_directory(str_dir, str_list) #loading required structures
                
                    structures2= align_structures(structures1) #aligning structures with the reference structure
                
                    avg_str= compute_average_structure(structures1, structures2) #calculating average structure
                
                    rmsds= align_avg(avg_str, structures1) #calculating distance from the average structure
                
                    rep_strs = min(rmsds, key=lambda k: rmsds[k])
                
                    for bp_IDs in rmsds:
                        df2.loc[(df2['bp_ID']== bp_IDs), "RMSD_to_avg"]= rmsds[bp_IDs] #check this line again
                
                    df2.loc[(df2['bp_ID']== rep_strs), "representative"]= 1
                    dfs.append(df2)
                
    df_r=pd.concat(dfs)
    df_r.index= np.arange(0, len(df_r))
    print (df_r.shape)
    return df_r        

#read the data
optparser = OptionParser()
(options, args) = optparser.parse_args()
csvfile = args[0]
D1= pd.read_csv(csvfile)

R= int(args[1])

#replace '.' in 'bp_ID' entries with '_'
D1['bp_ID'] = D1['bp_ID'].str.replace('.', '_')
D2= fix_PDB_ID(D1)

D3= rem_red(D2)

#store the finalized data as csv
if R==1:
    #for standard wobble
    D3.to_csv('../results/all_standard_wobble_redundancy_checked.csv', index= False)
elif R==2:
    #for shifted wobble
    D3.to_csv('../results/all_shifted_wobble_redundancy_checked.csv', index= False)
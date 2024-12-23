#importing required libraries
import os
import pandas as pd
from pymol import cmd, stored
from pymol import stored
from optparse import OptionParser
import numpy as np

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


def find_split_loc(x): 
    #this function will identify the location to split bp_res notation
    #bp_res example: A.G485-A.U499
    #sometimes the residue index can contain negative sign, which can add multiple '-' sign in the bp_res notation
    #those examples can be confusing for where to split
    ind=[]
    for k, l in enumerate(x):
        if l=='-':
            ind.append(k)
    if len(ind)>1:
        s=ind[1]
    else:
        s=ind[0]
    return s


def string_split(a):
    #function to extract chain ID (chain), residue name (resname), residue index (resnum), and icode
    #input (a) format: a.U1007 -->(chain ID).(residue ID)(residue index)
    b= list(a.split('.')[1])
    chain = list(a.split('.'))[0]
    resname = ''
    resnum = ''
    icode = '' #add an empty space here when this function is used to extract RSRZ values from RCSB dataset
    
    for i in range(len(b)):
        if b[i].isalpha():
            resname += b[i]
        else:
            break
        
    c= ''.join(b)
    d= c[len(resname):] #anything other than residue ID
    
    for i in range(len(d)):
        if d[i].isalpha():
            icode += d[i]
        else:
            resnum += d[i]
    
    if len(icode) ==1:
        pass
    else:
        icode = icode.replace(' ', '')
    #icode = icode.replace(' ', '')
    D={}
    D['chain'] = chain
    D['resname'] = resname
    D['resnum'] = resnum
    D['icode'] = icode
    
    return D


def tert_intr(df, n):
    hm_dr= os.getcwd()
    #df is the dataframe containing information for the base pair examples of the corresponding registers
    #dr is the directory where all the pdb/cif files are stored
    #n will indicate which register we are working with, for reg 1 --> n=1, for reg 2 --> n=2, for reg 3 --> n=3
    #df will be grouped by 'PDB_ID' and for each group corresponding 'bp_res' members will be stored as a list
    #adding empty columns which will be polulated with corresponding bond angles and distances
    df['G_N1_detail']=''
    df['G_N2_detail']=''
    df['G_N3_detail']=''
    df['G_O6_detail']=''
    df['G_N7_detail']=''
    df["G_O2'_detail"]=''
    df["G_O4'_detail"]=''
    df['U_O2_detail']=''
    df['U_N3_detail']=''
    df['U_O4_detail']=''
    df["U_O2'_detail"]=''
    df["U_O4'_detail"]=''

    df['G_N1_type']=''
    df['G_N2_type']=''
    df['G_N3_type']=''
    df['G_O6_type']=''
    df['G_N7_type']=''
    df["G_O2'_type"]=''
    df["G_O4'_type"]=''
    df['U_O2_type']=''
    df['U_N3_type']=''
    df['U_O4_type']=''
    df["U_O2'_type"]=''
    df["U_O4'_type"]=''

    df['combined_type']=''


    df1= df.groupby('PDB_ID')['bp_res'].apply(list).reset_index(name="list_bp_res")
    D={}
    for i, j in enumerate(df1['PDB_ID']):
        D[j]= df1['list_bp_res'][i]
    print (D)

    ##os.chdir(dr)
    work_stat= 0
    for p1 in D:
        #fname= p1+'.cif'
        #cmd.load(fname, fname) #corresponding pdb is loaded
        cmd.fetch(p1) #fetching the cif files directly from RCSB, no need to download or save in local directory
        cmd.remove('hydrogen') #all hydrogen atoms will be excluded
        #print (p1)
        for x1, x2 in enumerate(D[p1]):
            work_stat +=1
            print ('working here ------------------------------------------------->>'+ str(work_stat))
            print (p1)
            print (x2)
            a1= x2[:find_split_loc(x2)]
            b1= x2[find_split_loc(x2)+1:]

            if '^' in a1:
                a1 = a1.replace('^', '')
            else:
                pass
    
            if '^' in b1:
                b1 = b1.replace('^', '')
            else:
                pass
            
            a2 = string_split(a1) #chain, resname, resnum, icode for res1
            b2 = string_split(b1) #chain, resname, resnum, icode for res2

            a2['resnum']= a2['resnum']+a2['icode']
            b2['resnum']= b2['resnum']+b2['icode']

            a2['segi']= df[(df['PDB_ID']==p1)&(df['bp_res']==x2)]['seg_ID'].tolist()[0]
            b2['segi']= df[(df['PDB_ID']==p1)&(df['bp_res']==x2)]['seg_ID'].tolist()[0]

            if a2['resname']=='G':
                G_info= a2
                U_info= b2
            else:
                G_info= b2
                U_info= a2

            #getting the atom identifier numbers for the corresponding G and U atoms
            G_atoms= ["O4'", "O2'", 'N7', 'O6', 'N1', 'N2', 'N3'] #G atoms for which we want to identify atoms within 3.4 A
            G_ids= cmd.identify('(chain '+ G_info['chain'] + ' & '+ 'segi '+ G_info['segi'] + ' & ' +'resi '+ G_info['resnum']+  ") and (not name P+OP1+OP2+O5'+C5'+C4'+C3'+O3'+C2'+C1'+'N9'+'C8'+C5+C6+C2+C4)")
            G_ids.sort()
 
            U_atoms= ["O4'", "O2'", 'O2', 'N3', 'O4'] #U atoms for which we want to identify atoms within 3.4 A
            U_ids= cmd.identify('(chain '+ U_info['chain'] + ' & '+ 'segi '+ U_info['segi'] + ' & ' +'resi '+ U_info['resnum']+  ") and (not name P+OP1+OP2+O5'+C5'+C4'+C3'+O3'+C2'+C1'+N1+C2+C4+C5+C6)") 
            U_ids.sort()

            #now the calculation will move forward only if the residues are not missing any atoms
            if len(G_ids)==7 and len(U_ids)==5:
                G_ID= {G_atoms[g1]: str(G_ids[g1]) for g1 in range(len(G_atoms))}
                U_ID= {U_atoms[g1]: str(U_ids[g1]) for g1 in range(len(U_atoms))}
            
                #print (G_ID)
                #print (U_ID)


                bp_res_G_atoms= []
                bp_res_U_atoms= []
                bp_res_pairs= []

                if n==1:
                    #relevant atoms in G forming hydrogen bonds associated with the register of interest
                    bp_res_G_atoms.append(str(G_ID['O6'])+ '_'+ 'O6'+ '_'+ str(G_info["resnum"])+ '_'+ 'G')
                    bp_res_G_atoms.append(str(G_ID['N1'])+ '_'+ 'N1'+ '_'+ str(G_info["resnum"])+ '_'+ 'G')

                    #relevant atoms in U forming hydrogen bonds associated with the register of interest
                    bp_res_U_atoms.append(str(U_ID['N3'])+ '_'+ 'N3'+ '_'+ str(U_info["resnum"])+ '_'+ 'U')
                    bp_res_U_atoms.append(str(U_ID['O2'])+ '_'+ 'O2'+ '_'+ str(U_info["resnum"])+ '_'+ 'U')

                #elif n==2:
                #    #relevant atoms in G forming hydrogen bonds associated with the register of interest
                #    bp_res_G_atoms.append(str(G_ID['N1'])+ '_'+ 'N1'+ '_'+ str(G_info["resnum"])+ '_'+ 'G')
                #    bp_res_G_atoms.append(str(G_ID['N2'])+ '_'+ 'N2'+ '_'+ str(G_info["resnum"])+ '_'+ 'G')

                    #relevant atoms in U forming hydrogen bonds associated with the register of interest
                #    bp_res_U_atoms.append(str(U_ID['N3'])+ '_'+ 'N3'+ '_'+ str(U_info["resnum"])+ '_'+ 'U')
                #    bp_res_U_atoms.append(str(U_ID['O2'])+ '_'+ 'O2'+ '_'+ str(U_info["resnum"])+ '_'+ 'U')

                elif n==2:
                    #relevant atoms in G forming hydrogen bonds associated with the register of interest
                    bp_res_G_atoms.append(str(G_ID['N1'])+ '_'+ 'N1'+ '_'+ str(G_info["resnum"])+ '_'+ 'G')
                    bp_res_G_atoms.append(str(G_ID['N2'])+ '_'+ 'N2'+ '_'+ str(G_info["resnum"])+ '_'+ 'G')

                    #relevant atoms in U forming hydrogen bonds associated with the register of interest
                    bp_res_U_atoms.append(str(U_ID['O4'])+ '_'+ 'O4'+ '_'+ str(U_info["resnum"])+ '_'+ 'U')
                    bp_res_U_atoms.append(str(U_ID['N3'])+ '_'+ 'N3'+ '_'+ str(U_info["resnum"])+ '_'+ 'U')
                #print (bp_res_G_atoms)
                #print (bp_res_U_atoms)
                #bp_res_pairs.append()
                bp_res_pairs= [bp_res_G_atoms[0]+'-'+bp_res_U_atoms[0], 
                               bp_res_G_atoms[0]+'-'+bp_res_U_atoms[1], 
                               bp_res_G_atoms[1]+'-'+bp_res_U_atoms[0], 
                               bp_res_G_atoms[1]+'-'+bp_res_U_atoms[1],
                               bp_res_U_atoms[0]+'-'+bp_res_G_atoms[0], 
                               bp_res_U_atoms[0]+'-'+bp_res_G_atoms[1],
                               bp_res_U_atoms[1]+'-'+bp_res_G_atoms[0], 
                               bp_res_U_atoms[1]+'-'+bp_res_G_atoms[1]]
            
                #print (bp_res_pairs)
                total_atom_type=[]
                wres= {'G':[G_ID, G_info], 'U':[U_ID, U_info]}

                for res in wres:
                    for atms in wres[res][0]:
                        #print (atms)
                        atom_ID_det= wres[res][0][atms]+ '_'+ atms+ '_'+ str(wres[res][1]["resnum"])+ '_'+res
                        #below is the deistance wround hetero atoms to find tertiary/ non-WCF interactions
                        #future direction ====> add lines to ask user to define this distance while executing this script
                        all_surs= cmd.identify('(all within 3.4 of ID '+ str(wres[res][0][atms])+') and (not elem C) and (not elem H) and (not resi '+ str(wres[res][1]["resnum"])+')')

                        atom_detail=''
                        atom_type=''

                        for atm_in, atm_id in enumerate(all_surs):
                            stored.imp_sur4=[]
                            cmd.iterate("ID "+ str(atm_id), "stored.imp_sur4.append([name, resn, resi, chain, segi])")
                    
                            sur_atom_ID_det= str(atm_id)+ '_'+ stored.imp_sur4[0][0]+ '_'+ str(stored.imp_sur4[0][2])+ '_'+ stored.imp_sur4[0][1]
                            #print (sur_atom_ID_det)
                            interacting_pair= atom_ID_det+ '-'+ sur_atom_ID_det

                            if interacting_pair in bp_res_pairs:
                            #this interaction belong to the shifted wobble hydrogen bond
                                pass
                            else:
                                atom_name= stored.imp_sur4[0][1]
                                #print (atom_name)
                                if atom_name in residues:
                                    atom_type_in= residues[stored.imp_sur4[0][1]]['type']
                                else:
                                    atom_type_in= 'unknown'
                                if len(atom_type)==0:
                                    atom_type += atom_type_in
                                else:
                                    atom_type += '*'+atom_type_in
                        
                                if len(atom_detail) ==0:
                                    atom_detail += sur_atom_ID_det
                                else:
                                    atom_detail += '*'+ sur_atom_ID_det
                                total_atom_type.append(atom_type_in)
                        df.loc[((df['PDB_ID']==p1)&(df['bp_res']==x2)), res+"_"+ atms+ "_detail"] = atom_detail
                        df.loc[((df['PDB_ID']==p1)&(df['bp_res']==x2)), res+"_"+ atms+"_type"] = atom_type

                total_atom_type_str= '-'.join([a_type for a_type in list(set(total_atom_type)) if len(a_type)!=0])
                #print (total_atom_type_str)
                df.loc[((df['PDB_ID']==p1)&(df['bp_res']==x2)), "combined_type"] = total_atom_type_str
            
            else:
                df.loc[((df['PDB_ID']==p1)&(df['bp_res']==x2)), "combined_type"] = '####' #this value will help to identify observation with missing atoms




        cmd.delete(p1)
        cmd.delete('all')
        junk_cifs= [p1, p1.lower()]
        for jcifs in junk_cifs:
            jpaths= hm_dr+ '/'+ jcifs+ '.cif'
            if os.path.isfile(jpaths):
                os.remove(jpaths)
  
    junk_mcifs= ['UNK', 'K', 'MG', 'N']
    for jmcifs in junk_mcifs:
        jmpaths= hm_dr+ '/'+ jmcifs+ '.cif'
        if os.path.isfile(jmpaths):
            os.remove(jmpaths)
    df1= df[df['combined_type'] != '####'] #removing registers with missing atoms
    df1.index= np.arange(0, len(df1))

    return df1



residues = {
    "ARG": {
        "type": "amino_acid",
        "backbone_charge": "neutral",
        "side_chain_charge": "positive",
        "donors": {
            "name": ["N", "NE", "NH1", "NH2"],
            "trail": [["CA", "CB"], ["CZ", "NH1"], ["CZ", "NE"], ["CZ", "NE"]],
            "element": ["N", "N", "N", "N"],
            "hybridization": ["sp2", "sp2", "sp2", "sp2"],
            "rotatable": [False, False, False, False],
            "state": ["N", "N", "N", "N"]
        },
        "acceptors": {
            "name": ["O"],
            "trail": [["C", "CA"]],
            "element": ["O"],
            "hybridization": ["sp2"],
            "rotatable": [False],
            "state": ["N"]
        }
    },
    "LYS": {
        "type": "amino_acid",
        "backbone_charge": "neutral",
        "side_chain_charge": "positive",
        "donors": {
            "name": ["N", "NZ"],
            "trail": [["CA", "CB"], ["CE", "CD"]],
            "element": ["N", "N"],
            "hybridization": ["sp2", "sp3"],
            "rotatable": [False, True],
            "state": ["N", "N"]
        },
        "acceptors": {
            "name": ["O"],
            "trail": [["C", "CA"]],
            "element": ["O"],
            "hybridization": ["sp2"],
            "rotatable": [False],
            "state": ["N"]
        }
    },
    "HIS": {
        "type": "amino_acid",
        "backbone_charge": "neutral",
        "side_chain_charge": "positive",
        "donors": {
            "name": ["N", "ND1", "NE2"],
            "trail": [["CA", "CB"], ["CG", "CD2"], ["CD2", "CG"]],
            "element": ["N", "N", "N"],
            "hybridization": ["sp2", "sp2", "sp2"],
            "rotatable": [False, False, False],
            "state": ["N", "N", "N"]
        },
        "acceptors": {
            "name": ["O", "ND1", "NE2"],
            "trail": [["C", "CA"], ["CG", "CD2"], ["CD2", "CG"]],
            "element": ["O", "N", "N"],
            "hybridization": ["sp2", "sp2", "sp2"],
            "rotatable": [False, False, False],
            "state": ["N", "N", "N"]
        }
    },
    "TRP": {
        "type": "amino_acid",
        "backbone_charge": "neutral",
        "side_chain_charge": "neutral",
        "donors": {
            "name": ["N", "NE1"],
            "trail": [["CA", "CB"], ["CD1", "CG"]],
            "element": ["N", "N"],
            "hybridization": ["sp2", "sp2"],
            "rotatable": [False, False],
            "state": ["N", "N"]
        },
        "acceptors": {
            "name": ["O"],
            "trail": [["C", "CA"]],
            "element": ["O"],
            "hybridization": ["sp2"],
            "rotatable": [False],
            "state": ["N"]
        }
    },
    "ASP": {
        "type": "amino_acid",
        "backbone_charge": "neutral",
        "side_chain_charge": "negative",
        "donors": {
            "name": ["N"],
            "trail": [["CA", "CB"]],
            "element": ["N"],
            "hybridization": ["sp2"],
            "rotatable": [False],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["O", "OD1", "OD2"],
            "trail": [["C", "CA"], ["CG", "CB"], ["CG", "CB"]],
            "element": ["O", "O", "O"],
            "hybridization": ["sp2", "sp2", "sp2"],
            "rotatable": [False, False, False],
            "state": ["N", "N", "N"]
        }
    },
    "GLU": {
        "type": "amino_acid",
        "backbone_charge": "neutral",
        "side_chain_charge": "negative",
        "donors": {
            "name": ["N"],
            "trail": [["CA", "CB"]],
            "element": ["N"],
            "hybridization": ["sp2"],
            "rotatable": [False],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["O", "OE1", "OE2"],
            "trail": [["C", "CA"], ["CD", "CG"], ["CD", "CG"]],
            "element": ["O", "O", "O"],
            "hybridization": ["sp2", "sp2", "sp2"],
            "rotatable": [False, False, False],
            "state": ["N", "N", "N"]
        }
    },
    "SER": {
        "type": "amino_acid",
        "backbone_charge": "neutral",
        "side_chain_charge": "neutral",
        "donors": {
            "name": ["N", "OG"],
            "trail": [["CA", "CB"], ["CB", "CA"]],
            "element": ["N", "O"],
            "hybridization": ["sp2", "sp3"],
            "rotatable": [False, True],
            "state": ["N", "N"]
        },
        "acceptors": {
            "name": ["O", "OG"],
            "trail": [["C", "CA"], ["CB", "CA"]],
            "element": ["O", "O"],
            "hybridization": ["sp2", "sp3"],
            "rotatable": [False, True],
            "state": ["N", "N"]
        },
    },
    "CYS": {
        "type": "amino_acid",
        "backbone_charge": "neutral",
        "side_chain_charge": "neutral",
        "donors": {
            "name": ["N", "SG"],
            "trail": [["CA", "CB"], ["CB", "CA"]],
            "element": ["N", "S"],
            "hybridization": ["sp2", "sp3"],
            "rotatable": [False, True],
            "state": ["N", "N"]
        },
        "acceptors": {
            "name": ["O", "SG"],
            "trail": [["C", "CA"], ["CB", "CA"]],
            "element": ["O", "S"],
            "hybridization": ["sp2", "sp3"],
            "rotatable": [False, True],
            "state": ["N", "N"]
        },
    },
    "MET": {
        "type": "amino_acid",
        "backbone_charge": "neutral",
        "side_chain_charge": "neutral",
        "donors": {
            "name": ["N"],
            "trail": [["CA", "CB"]],
            "element": ["N"],
            "hybridization": ["sp2"],
            "rotatable": [False],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["O", "SD"],
            "trail": [["C", "CA"], ["CG", "CB"]],
            "element": ["O", "S"],
            "hybridization": ["sp2", "sp3"],
            "rotatable": [False, False],
            "state": ["N", "N"]
        },
    },
    "THR": {
        "type": "amino_acid",
        "backbone_charge": "neutral",
        "side_chain_charge": "neutral",
        "donors": {
            "name": ["N", "OG1"],
            "trail": [["CA", "CB"], ["CB", "CA"]],
            "element": ["N", "O"],
            "hybridization": ["sp2", "sp3"],
            "rotatable": [False, True],
            "state": ["N", "N"]
        },
        "acceptors": {
            "name": ["O", "OG1"],
            "trail": [["C", "CA"], ["CB", "CA"]],
            "element": ["O", "O"],
            "hybridization": ["sp2", "sp3"],
            "rotatable": [False, True],
            "state": ["N", "N"]
        },
    },
    "ASN": {
        "type": "amino_acid",
        "backbone_charge": "neutral",
        "side_chain_charge": "negative",
        "donors": {
            "name": ["N", "ND2"],
            "trail": [["CA", "CB"], ["CG", "CB"]],
            "element": ["N", "N"],
            "hybridization": ["sp2", "sp2"],
            "rotatable": [False, False],
            "state": ["N", "N"]
        },
        "acceptors": {
            "name": ["O", "OD1"],
            "trail": [["C", "CA"], ["CG", "CB"]],
            "element": ["O", "O"],
            "hybridization": ["sp2", "sp2"],
            "rotatable": [False, False],
            "state": ["N", "N"]
        },
    },
    "GLN": {
        "type": "amino_acid",
        "backbone_charge": "neutral",
        "side_chain_charge": "negative",
        "donors": {
            "name": ["N", "NE2"],
            "trail": [["CA", "CB"], ["CD", "CG"]],
            "element": ["N", "N"],
            "hybridization": ["sp2", "sp2"],
            "rotatable": [False, False],
            "state": ["N", "N"]
        },
        "acceptors": {
            "name": ["O", "OE1"],
            "trail": [["C", "CA"], ["CD", "CG"]],
            "element": ["O", "O"],
            "hybridization": ["sp2", "sp2"],
            "rotatable": [False, False],
            "state": ["N", "N"]
        },
    },
    "ALA": {
        "type": "amino_acid",
        "backbone_charge": "neutral",
        "side_chain_charge": "neutral",
        "donors": {
            "name": ["N"],
            "trail": [["CA", "CB"]],
            "element": ["N"],
            "hybridization": ["sp2"],
            "rotatable": [False],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["O"],
            "trail": [["C", "CA"]],
            "element": ["O"],
            "hybridization": ["sp2"],
            "rotatable": [False],
            "state": ["N"]
        },
    },
    "GLY": {
        "type": "amino_acid",
        "backbone_charge": "neutral",
        "side_chain_charge": "neutral",
        "donors": {
            "name": ["N"],
            "trail": [["CA", "C"]],
            "element": ["N"],
            "hybridization": ["sp2"],
            "rotatable": [False],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["O"],
            "trail": [["C", "CA"]],
            "element": ["O"],
            "hybridization": ["sp2"],
            "rotatable": [False],
            "state": ["N"]
        },
    },
    "ILE": {
        "type": "amino_acid",
        "backbone_charge": "neutral",
        "side_chain_charge": "neutral",
        "donors": {
            "name": ["N"],
            "trail": [["CA", "CB"]],
            "element": ["N"],
            "hybridization": ["sp2"],
            "rotatable": [False],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["O"],
            "trail": [["C", "CA"]],
            "element": ["O"],
            "hybridization": ["sp2"],
            "rotatable": [False],
            "state": ["N"]
        },
    },
    "LEU": {
        "type": "amino_acid",
        "backbone_charge": "neutral",
        "side_chain_charge": "neutral",
        "donors": {
            "name": ["N"],
            "trail": [["CA", "CB"]],
            "element": ["N"],
            "hybridization": ["sp2"],
            "rotatable": [False],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["O"],
            "trail": [["C", "CA"]],
            "element": ["O"],
            "hybridization": ["sp2"],
            "rotatable": [False],
            "state": ["N"]
        },
    },
    "PRO": {
        "type": "amino_acid",
        "backbone_charge": "neutral",
        "side_chain_charge": "neutral",
        "donors": {
            "name": [],
            "trail": [],
            "element": [],
            "hybridization": [],
            "rotatable": [],
            "state": []
        },
        "acceptors": {
            "name": ["O"],
            "trail": [["C", "CA"]],
            "element": ["O"],
            "hybridization": ["sp2"],
            "rotatable": [False],
            "state": ["N"]
        },
    },
    "PHE": {
        "type": "amino_acid",
        "backbone_charge": "neutral",
        "side_chain_charge": "neutral",
        "donors": {
            "name": ["N"],
            "trail": [["CA", "CB"]],
            "element": ["N"],
            "hybridization": ["sp2"],
            "rotatable": [False],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["O"],
            "trail": [["C", "CA"]],
            "element": ["O"],
            "hybridization": ["sp2"],
            "rotatable": [False],
            "state": ["N"]
        },
    },
    "TYR": {
        "type": "amino_acid",
        "backbone_charge": "neutral",
        "side_chain_charge": "neutral",
        "donors": {
            "name": ["N", "OH"],
            "trail": [["CA", "CB"], ["CZ", "CE1"]],
            "element": ["N", "O"],
            "hybridization": ["sp2", "sp3"],
            "rotatable": [False, True],
            "state": ["N", "N"]
        },
        "acceptors": {
            "name": ["O", "OH"],
            "trail": [["C", "CA"], ["CZ", "CE1"]],
            "element": ["O", "O"],
            "hybridization": ["sp2", "sp3"],
            "rotatable": [False, True],
            "state": ["N", "N"]
        },
    },
    "VAL": {
        "type": "amino_acid",
        "backbone_charge": "neutral",
        "side_chain_charge": "neutral",
        "donors": {
            "name": ["N"],
            "trail": [["CA", "CB"]],
            "element": ["N"],
            "hybridization": ["sp2"],
            "rotatable": [False],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["O"],
            "trail": [["C", "CA"]],
            "element": ["O"],
            "hybridization": ["sp2"],
            "rotatable": [False],
            "state": ["N"]
        },
    },
    "G": {
        "type": "RNA",
        "backbone_charge": "negative",
        "side_chain_charge": "neutral",
        "donors": {
            "name": ["O2'", "N1", "N2", "O6"],
            "trail": [["C2'", "C1'"], ["C2", "N3"], ["C2", "N3"], ["C6", "C5"]],
            "element": ["O", "N", "N", "O"],
            "hybridization": ["sp3", "sp2", "sp2", "sp2"],
            "rotatable": [True, False, False, False],
            "state": ["N", "N", "N", "T"]
        },
        "acceptors": {
            "name": ["OP1", "OP2", "O5'", "O4'", "O3'", "O2'", "N3", "O6", "N7", "N1"],
            "trail": [["P", "O5'"], ["P", "O5'"], ["C5'", "C4'"], ["C4'", "C3'"], ["C3'", "C2'"], ["C2'", "C1'"],
                      ["C4", "C5"], ["C6", "N1"], ["C8", "N9"], ["C2", "N3"]],
            "element": ["O", "O", "O", "O", "O", "O", "N", "O", "N", "N"],
            "hybridization": ["sp2", "sp2", "sp3", "sp3", "sp3", "sp3", "sp2", "sp2", "sp2", "sp2"],
            "rotatable": [False, False, False, False, False, True, False, False, False, False],
            "state": ["N", "N", "N", "N", "N", "N", "N", "N", "N", "T"]
        }
    },
    "A": {
        "type": "RNA",
        "backbone_charge": "negative",
        "side_chain_charge": "neutral",
        "donors": {
            "name": ["O2'", "N6"],
            "trail": [["C2'", "C1'"], ["C6", "C5"]],
            "element": ["O", "N"],
            "hybridization": ["sp3", "sp2"],
            "rotatable": [True, False],
            "state": ["N", "N"]
        },
        "acceptors": {
            "name": ["OP1", "OP2", "O5'", "O4'", "O3'", "O2'", "N1", "N3", "N7"],
            "trail": [["P", "O5'"], ["P", "O5'"], ["C5'", "C4'"], ["C4'", "C3'"], ["C3'", "C2'"], ["C2'", "C1'"],
                      ["C2", "N3"], ["C2", "N1"], ["C5", "C4"]],
            "element": ["O", "O", "O", "O", "O", "O", "N", "N", "N"],
            "hybridization": ["sp2", "sp2", "sp3", "sp3", "sp3", "sp3", "sp2", "sp2", "sp2"],
            "rotatable": [False, False, False, False, False, True, False, False, False],
            "state": ["N", "N", "N", "N", "N", "N", "N", "N", "N"]
        }
    },
    "C": {
        "type": "RNA",
        "backbone_charge": "negative",
        "side_chain_charge": "neutral",
        "donors": {
            "name": ["O2'", "N4"],
            "trail": [["C2'", "C1'"], ["C4", "C5"]],
            "element": ["O", "N"],
            "hybridization": ["sp3", "sp2"],
            "rotatable": [True, False],
            "state": ["N", "N"]
        },
        "acceptors": {
            "name": ["OP1", "OP2", "O5'", "O4'", "O3'", "O2'", "N3", "O2"],
            "trail": [["P", "O5'"], ["P", "O5'"], ["C5'", "C4'"], ["C4'", "C3'"], ["C3'", "C2'"], ["C2'", "C1'"],
                      ["C2", "N1"], ["C2", "N1"]],
            "element": ["O", "O", "O", "O", "O", "O", "N", "O"],
            "hybridization": ["sp2", "sp2", "sp3", "sp3", "sp3", "sp3", "sp2", "sp2"],
            "rotatable": [False, False, False, False, False, True, False, False],
            "state": ["N", "N", "N", "N", "N", "N", "N", "N"]
        }
    },
    "U": {
        "type": "RNA",
        "backbone_charge": "negative",
        "side_chain_charge": "neutral",
        "donors": {
            "name": ["O2'", "N3", "O2", "O4"],
            "trail": [["C2'", "C1'"], ["C2", "N1"], ["C2", "N1"], ["C4", "N3"]],
            "element": ["O", "N", "O", "O"],
            "hybridization": ["sp3", "sp2", "sp2", "sp2"],
            "rotatable": [True, False, False, False],
            "state": ["N", "N", "T", "T"]
        },
        "acceptors": {
            "name": ["OP1", "OP2", "O5'", "O4'", "O3'", "O2'", "O4", "O2", "N3"],
            "trail": [["P", "O5'"], ["P", "O5'"], ["C5'", "C4'"], ["C4'", "C3'"],
                      ["C3'", "C2'"], ["C2'", "C1'"], ["C4", "C5"], ["C2", "N3"], ["C2", "N1"]],
            "element": ["O", "O", "O", "O", "O", "O", "O", "O", "N"],
            "hybridization": ["sp2", "sp2", "sp3", "sp3", "sp3", "sp3", "sp2", "sp2", "sp2"],
            "rotatable": [False, False, False, False, False, True, False, False, False],
            "state": ["N", "N", "N", "N", "N", "N", "N", "N", "T"]
        }
    },
    "DG": {
        "type": "DNA",
        "backbone_charge": "negative",
        "side_chain_charge": "neutral",
        "donors": {
            "name": ["N1", "N2", "O6"],
            "trail": [["C2", "N3"], ["C2", "N3"], ["C6", "C5"]],
            "element": ["N", "N", "O"],
            "hybridization": ["sp2", "sp2", "sp2"],
            "rotatable": [False, False, False],
            "state": ["N", "N", "T"]
        },
        "acceptors": {
            "name": ["OP1", "OP2", "O5'", "O4'", "O3'", "N3", "O6", "N7", "N1"],
            "trail": [["P", "O5'"], ["P", "O5'"], ["C5'", "C4'"], ["C4'", "C3'"], ["C3'", "C2'"], ["C4", "C5"],
                      ["C6", "N1"], ["C8", "N9"], ["C2", "N3"]],
            "element": ["O", "O", "O", "O", "O", "N", "O", "N", "N"],
            "hybridization": ["sp2", "sp2", "sp3", "sp3", "sp3", "sp2", "sp2", "sp2", "sp2"],
            "rotatable": [False, False, False, False, False, False, False, False, False],
            "state": ["N", "N", "N", "N", "N", "N", "N", "N", "T"]
        }
    },
    "DA": {
        "type": "DNA",
        "backbone_charge": "negative",
        "side_chain_charge": "neutral",
        "donors": {
            "name": ["N6"],
            "trail": [["C6", "C5"]],
            "element": ["N"],
            "hybridization": ["sp2"],
            "rotatable": [False],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["OP1", "OP2", "O5'", "O4'", "O3'", "N1", "N3", "N7"],
            "trail": [["P", "O5'"], ["P", "O5'"], ["C5'", "C4'"], ["C4'", "C3'"], ["C3'", "C2'"], ["C2", "N3"],
                      ["C2", "N1"], ["C5", "C4"]],
            "element": ["O", "O", "O", "O", "O", "N", "N", "N"],
            "hybridization": ["sp2", "sp2", "sp3", "sp3", "sp3", "sp2", "sp2", "sp2"],
            "rotatable": [False, False, False, False, False, False, False, False],
            "state": ["N", "N", "N", "N", "N", "N", "N", "N"]
        }
    },
    "DC": {
        "type": "DNA",
        "backbone_charge": "negative",
        "side_chain_charge": "neutral",
        "donors": {
            "name": ["N4"],
            "trail": [["C4", "C5"]],
            "element": ["N"],
            "hybridization": ["sp2"],
            "rotatable": [False],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["OP1", "OP2", "O5'", "O4'", "O3'", "N3", "O2"],
            "trail": [["P", "O5'"], ["P", "O5'"], ["C5'", "C4'"], ["C4'", "C3'"], ["C3'", "C2'"], ["C2", "N1"],
                      ["C2", "N1"]],
            "element": ["O", "O", "O", "O", "O", "N", "O"],
            "hybridization": ["sp2", "sp2", "sp3", "sp3", "sp3", "sp2", "sp2"],
            "rotatable": [False, False, False, False, False, False, False],
            "state": ["N", "N", "N", "N", "N", "N", "N"]
        }
    },
    "DT": {
        "type": "DNA",
        "backbone_charge": "negative",
        "side_chain_charge": "neutral",
        "donors": {
            "name": ["N3", "O2", "O4"],
            "trail": [["C2", "N1"], ["C2", "N1"], ["C4", "N3"]],
            "element": ["N", "O", "O"],
            "hybridization": ["sp2", "sp2", "sp2"],
            "rotatable": [False, False, False],
            "state": ["N", "T", "T"]
        },
        "acceptors": {
            "name": ["OP1", "OP2", "O5'", "O4'", "O3'", "O4", "O2", "N3"],
            "trail": [["P", "O5'"], ["P", "O5'"], ["C5'", "C4'"], ["C4'", "C3'"],
                      ["C3'", "C2'"], ["C4", "C5"], ["C2", "N3"], ["C2, N1"]],
            "element": ["O", "O", "O", "O", "O", "O", "O", "N"],
            "hybridization": ["sp2", "sp2", "sp3", "sp3", "sp3", "sp2", "sp2", "sp2"],
            "rotatable": [False, False, False, False, False, False, False, False],
            "state": ["N", "N", "N", "N", "N", "N", "N", "T"]
        }
    },
    "HOH": {
        "type": "solvent",
        "donors": {
            "name": ["O"],
            "element": ["O"],
            "hybridization": ["sp3"],
            "rotatable": [True],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["O"],
            "element": ["O"],
            "hybridization": ["sp3"],
            "rotatable": [True],
            "state": ["N"]
        }
    },
    "NA": {
        "type": "ion",
        "charge": 1
    },
    "K": {
        "type": "ion",
        "charge": 1
    },
    "TL": {
        "type": "ion",
        "charge": 1
    },
    "CU1": {
        "type": "ion",
        "charge": 1
    },
    "AG": {
        "type": "ion",
        "charge": 1
    },
    "RH": {
        "type": "ion",
        "charge": 1
    },
    "AU": {
        "type": "ion",
        "charge": 1
    },
    "MG": {
        "type": "ion",
        "charge": 2
    },
    "MN": {
        "type": "ion",
        "charge": 2
    },
    "0BE": {
        "type": "ion",
        "charge": 2
    },
    "BA": {
        "type": "ion",
        "charge": 2
    },
    "CA": {
        "type": "ion",
        "charge": 2
    },
    "SR": {
        "type": "ion",
        "charge": 2
    },
    "CO": {
        "type": "ion",
        "charge": 2
    },
    "CU": {
        "type": "ion",
        "charge": 2
    },
    "FE2": {
        "type": "ion",
        "charge": 2
    },
    "NI": {
        "type": "ion",
        "charge": 2
    },
    "ZN": {
        "type": "ion",
        "charge": 2
    },
    "CD": {
        "type": "ion",
        "charge": 2
    },
    "PD": {
        "type": "ion",
        "charge": 2
    },
    "Y1": {
        "type": "ion",
        "charge": 2
    },
    "HG": {
        "type": "ion",
        "charge": 2
    },
    "PT": {
        "type": "ion",
        "charge": 2
    },
    "PB": {
        "type": "ion",
        "charge": 2
    },
    "3CO": {
        "type": "ion",
        "charge": 3
    },
    "3NI": {
        "type": "ion",
        "charge": 3
    },
    "CR": {
        "type": "ion",
        "charge": 3
    },
    "CU3": {
        "type": "ion",
        "charge": 3
    },
    "FE": {
        "type": "ion",
        "charge": 3
    },
    "MN3": {
        "type": "ion",
        "charge": 3
    },
    "V": {
        "type": "ion",
        "charge": 3
    },
    "RH3": {
        "type": "ion",
        "charge": 3
    },
    "RU": {
        "type": "ion",
        "charge": 3
    },
    "YT3": {
        "type": "ion",
        "charge": 3
    },
    "AU3": {
        "type": "ion",
        "charge": 3
    },
    "IR3": {
        "type": "ion",
        "charge": 3
    },
    "LA": {
        "type": "ion",
        "charge": 3
    },
    "OS": {
        "type": "ion",
        "charge": 3
    },
    "SB": {
        "type": "ion",
        "charge": 3
    },
    "AL": {
        "type": "ion",
        "charge": 3
    },
    "BS3": {
        "type": "ion",
        "charge": 3
    },
    "IN": {
        "type": "ion",
        "charge": 3
    },
    "4TI": {
        "type": "ion",
        "charge": 4
    },
    "4MO": {
        "type": "ion",
        "charge": 4
    },
    "ZR": {
        "type": "ion",
        "charge": 4
    },
    "IR": {
        "type": "ion",
        "charge": 4
    },
    "OS4": {
        "type": "ion",
        "charge": 4
    },
    "PT4": {
        "type": "ion",
        "charge": 4
    },
    "6MO": {
        "type": "ion",
        "charge": 6
    },
    "W": {
        "type": "ion",
        "charge": 6
    },
    "NCO": {
        "type": "ion",
        "charge": 3,
        "donors": {
            "name": ["N1", "N2", "N3", "N4", "N5", "N6"],
            "trail": [["CO", "N3"], ["CO", "N4"], ["CO", "N1"], ["CO", "N2"], ["CO", "N6"], ["CO", "N5"]],
            "element": ["N", "N", "N", "N", "N", "N"],
            "hybridization": ["sp3", "sp3", "sp3", "sp3", "sp3", "sp3"],
            "rotatable": [True, True, True, True, True, True],
            "state": ["N", "N", "N", "N", "N", "N"]
        }
    }
}


optparser = OptionParser()
(options, args) = optparser.parse_args()
csvfile = args[0]
D1= pd.read_csv(csvfile, keep_default_na=False, dtype={'seg_ID': str, 'res_ID_res1': str, 'seg_ID': str, 'res_ID_res2': str})
 
R= int(args[1]) #R is indicating the register we are analyzing here
D2= tert_intr(D1, R)

if R==1:
    #for standard wobble
    D2.to_csv('../results/all_standard_wobble_with_non_WCF_interactions.csv', index= False)
elif R==2:
    #for shifted wobble
    D2.to_csv('../results/all_shifted_wobble_with_non_WCF_interaction.csv', index= False)


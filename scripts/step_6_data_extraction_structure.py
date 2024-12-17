#importing required libraries
import os
import pandas as pd
from pymol import cmd
from pymol import stored
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

def temp_factor(p, c, s, sp, dr):
    #p= PDB_ID
    #c= chain ID
    #s= segment ID
    #sp= can be '0' or '1' --> 0= no atom will be excluded, 1= exclude only sugar atoms, 2= exclude only phosphate atoms, 3= exclude both sugar and phosphate atoms
    #dr= direcory where all the structure files are stored
    
    #notes on Jan 17, 2024: mother function is already changing the directory, therefore, no need to do this here
    #hm_dr= os.getcwd() #this is the directory where this script is stored
    
    #changing the directory to the location where all structure files are stored
    #notes on Jan 17, 2024: mother function is already changing the directory, therefore, no need to do this here
    #os.chdir(dr)
    
    #loading the corresponding structure files
    #notes on Nov 15, 2023: will be loaded in the mother function (dis_ang), therefore, no need to load again
    #fname= p+'.cif'
    #cmd.load(fname, fname)

    #storeing all residues of the chain of interest
    stored.seq={} #in this dictionary key will be residue index and value will residue ID/name
    cmd.iterate('chain '+ c + ' and segi '+ s, 'stored.seq[resi]= resn')
    
    #storing average temp factor for each residues
    atom_bs={} #in this dictionary key will be residue index and value will be average temp factor of the atoms of corresponding residue

    for i in stored.seq:
        stored.bs={} #in this dictionary, key will be the atom identifier number, and 
        print (i)
        print (stored.seq[i])

        if sp==0:
            x= (cmd.identify('(chain '+ c + ' & '+ 'segi '+ s + ' & ' +'resi '+ str(i)+ ' and resn '+ stored.seq[i]+ ")")) #neither sugar nor phosphate atoms are excluded
        elif sp == 1:
            x= (cmd.identify('(chain '+ c + ' & '+ 'segi '+ s + ' & ' +'resi '+ str(i)+ ' and resn '+ stored.seq[i]+ ") and (not name O5'+C5'+C4'+O4'+C3'+O3'+C2'+O2'+C1')")) #only sugar atoms are excluded
        elif sp == 2:
            x= (cmd.identify('(chain '+ c + ' & '+ 'segi '+ s + ' & ' +'resi '+ str(i)+ ' and resn '+ stored.seq[i]+ ") and (not name P+OP1+OP2)")) #only phosphate atoms are excluded
        elif sp == 3:
            x= (cmd.identify('(chain '+ c + ' & '+ 'segi '+ s + ' & ' +'resi '+ str(i)+ ' and resn '+ stored.seq[i]+ ") and (not name P+OP1+OP2+O5'+C5'+C4'+O4'+C3'+O3'+C2'+O2'+C1')")) #both sugar and phosphate atoms are excluded
        
        y= sorted(x) #sorting the indices of the atoms from the corresponding residues 

        for j, k in enumerate(y):
            cmd.iterate('chain '+ c + ' and segi '+ s + ' and resi '+ str(i) + ' and resn '+ stored.seq[i] + 'and ID ' +str(k), 'stored.bs[ID]=b')
        #print (stored.bs)
        bs=[]
        for l in stored.bs:
            bs.append(stored.bs[l])
        print (bs)
        if len(bs)==0:
            print ('NO RESIDUES ARE FOUND')
            pass            
        else:
            avg_b= (sum(bs)/len(bs))
            print (avg_b)
            atom_bs[i]= avg_b
        #print (y)

    #notes on Nov 15, 2023: more calculation will be done in the mother function (dis_ang), therefore, no need to close in this function     
    #cmd.delete(fname) #closing/ deleting the opend structure file
    
    #notes on Jan 17, 2024: mother function is already changing the directory, therefore, no need to do this here
    #os.chdir(hm_dr) #going back to the home directory
    return atom_bs

def get_info_structure(df, n):
    #updated on February 8, 2024
    #finction to calculate required bond distances and dihedral angles and temperature factors
    #input will be a csv file with column 'bp_res', and the directory containing all the cif files
    #data format in 'bp_res' column: A.G2063-A.U2439
    #sp will indicate which atoms to consider for temperature factor calculation

    hm_dr= os.getcwd() #home directory where this script and the df is stored
    
    #adding empty columns which will be polulated with corresponding bond angles and distances
    #notes on June 13, 2024
    #add 9 distance columns
    WCF_G= ['N1', 'N2', 'O6']
    WCF_U= ['O2', 'N3', 'O4']
    #new lines from June 13, 2024 ends here

    ##df['b0_dis']= ''
    ##df['b0_dih']= ''
    ##df['b1_dis']= ''
    ##df['b1_dih']= ''
    ##df['b2_dis']= ''
    ##df['b2_dih']= '' #new directionality: from U to G
    ##df["C1'_dis"]= ''
    ##df["C1'_dih"]= ''

    #adding columns for the length of chain in structure
    ##df.insert (df.columns.get_loc('Chain_length_reference')+1, 'Chain_length_structure', '')
    
    #adding columns for temperature factors of G, U, and the whole molecule
    df['b_resG'] = ''
    df['b_resU'] = ''
    df['b_comb'] = ''
    ####df['b_all'] = ''
    
    #df will be grouped by 'PDB_ID' and for each group corresponding 'bp_res' members will be stored as a list
    df1= df.groupby('PDB_ID')['bp_res'].apply(list).reset_index(name="list_bp_res")
    D={}
    for i, j in enumerate(df1['PDB_ID']):
        D[j]= df1['list_bp_res'][i]
    print (D)
    
    status_count=0
    for p1 in D:
        #fname= p1+'.cif'
        #cmd.load(fname, fname) #corresponding pdb is loaded
        cmd.fetch(p1) #fetching the cif files directly from RCSB, no need to download or save in local directory
        cmd.remove('hydrogen') #all hydrogen atoms will be excluded
        print (p1)
        for x1, x2 in enumerate(D[p1]):
            status_count +=1
            print ('WORKING HERE --------------------------------------------------------------------------------------------------------------------- '+ str(x2))
            print (status_count)
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
            #extracting temperature factor 
            #storing average temp factor for each residues
            #in the all_bs dictionary, key will be the residue index and value will be the corresponding average temperature factor
            #in all_bs, to calculate the temperature factors, both sugar and phosphate residues were excluded 
                
            ####all_bs= temp_factor(p1, a2['chain'], a2['segi'], 3, dr) #chain ID and segment ID is same for both G and U. The input parameter '3' (sp as function input) can be added here or as function input
            ####G_b= all_bs[G_info['resnum']]
            ####U_b= all_bs[U_info['resnum']]
            ####b_comb= (G_b+U_b)/2 #avg temperature factor of G and U
            ####b_all= (sum(list(all_bs[rs] for rs in all_bs))/len(all_bs))
                
            #getting the atom identifier numbers for the corresponding G and U atoms
            G_atoms= ["C1'", 'N9', 'C8', 'N7', 'C5', 'C6', 'O6', 'N1', 'C2', 'N2', 'N3', 'C4'] #updated the order of atoms on Feb 8, 2023 
            G_ids= cmd.identify('(chain '+ G_info['chain'] + ' & '+ 'segi '+ G_info['segi'] + ' & ' +'resi '+ G_info['resnum']+  ") and (not name P+OP1+OP2+O5'+C5'+C4'+O4'+C3'+O3'+C2'+O2')")
            G_ids.sort()

            U_atoms= ["C1'", 'N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6'] #updated the order of atoms on Feb 8, 2023
            U_ids= cmd.identify('(chain '+ U_info['chain'] + ' & '+ 'segi '+ U_info['segi'] + ' & ' +'resi '+ U_info['resnum']+  ") and (not name P+OP1+OP2+O5'+C5'+C4'+O4'+C3'+O3'+C2'+O2')") 
            U_ids.sort()
            
            #now the calculation will move forward only if the residues are not missing any atoms
            if len(G_ids)==12 and len(U_ids)==9:
                G_ID= {G_atoms[g1]: str(G_ids[g1]) for g1 in range(len(G_atoms))}
                U_ID= {U_atoms[g1]: str(U_ids[g1]) for g1 in range(len(U_atoms))}

                #calculating distances and dihedral angles common for different wobbles 
                #b0= O6-O4
                #b0_dis= (cmd.get_distance(atom1= 'ID ' +G_ID['O6'], atom2= 'ID '+U_ID['O4'])) #distance between G.O6 and U.O4
                O6_O4_dih= (cmd.get_dihedral(atom1= 'ID ' +G_ID['C6'], atom2= 'ID '+G_ID['O6'], atom3= 'ID '+U_ID['O4'], atom4= 'ID '+U_ID['C4'])) #dihedral angle: G.C6-G.O6-U.O4-U.C4
                df.loc[((df['PDB_ID']==p1)&(df['bp_res']==x2)), "O6_O4_dih"] = O6_O4_dih
            
                #followings are distances and dihedral angles for C1' sugar atoms for G and U
                C1_dis= (cmd.get_distance(atom1= 'ID ' +G_ID["C1'"], atom2= 'ID '+U_ID["C1'"])) #distance between G.C1' and U.C1'
                df.loc[((df['PDB_ID']==p1)&(df['bp_res']==x2)), "C1_dis"] = C1_dis
                C1_dih = (cmd.get_dihedral(atom1= 'ID ' +G_ID['N9'], atom2= 'ID '+G_ID["C1'"], atom3= 'ID '+U_ID["C1'"], atom4= 'ID '+U_ID['N1'])) #dihedral angle: G.N9-G.C1'-U.C1'-U.N1
                df.loc[((df['PDB_ID']==p1)&(df['bp_res']==x2)), "C1_dih"] = C1_dih
            
                #notes on June 13, 2024
                #followings are distances and dihedral angles different for different rare wobbles
                #calculate all-against-all WCF face atom distance for G and U
                distances= {}
                for atomi, atomj in enumerate(WCF_G):
                    for atomk, atoml in enumerate(WCF_U):
                        #calculate the dihedral angles here
                        dis_name= atomj+ '_'+ atoml+ '_dis'
                        df.loc[((df['PDB_ID']==p1)&(df['bp_res']==x2)), dis_name]= (cmd.get_distance(atom1= 'ID ' +G_ID[atomj], atom2= 'ID '+U_ID[atoml]))
                        #distances[dis_name]= (cmd.get_distance(atom1= 'ID ' +G_ID[atomj], atom2= 'ID '+U_ID[atoml]))
                #new lines from June 13, 2024 ends here


                if n==1: #standard G•U wobble
                    #b1= O6-N3
                    #b1_dis= (cmd.get_distance(atom1= 'ID ' +G_ID['O6'], atom2= 'ID '+U_ID['N3'])) #distance between G.O6 and U.N3
                    O6_N3_dih= (cmd.get_dihedral(atom1= 'ID ' +G_ID['C6'], atom2= 'ID ' +G_ID['O6'], atom3= 'ID '+U_ID['N3'], atom4= 'ID '+U_ID['C2'])) #dihedral angle: G.C6-G.O6-U.N3-U.C2
                    df.loc[((df['PDB_ID']==p1)&(df['bp_res']==x2)), "O6_N3_dih"]= O6_N3_dih
                    #b2= N1-O2
                    #b2_dis= (cmd.get_distance(atom1= 'ID ' +G_ID['N1'], atom2= 'ID '+U_ID['O2'])) #distance between G.N1 and U.O2
                    N1_O2_dih= (cmd.get_dihedral(atom1= 'ID ' +G_ID['C6'], atom2= 'ID '+G_ID['N1'], atom3= 'ID '+U_ID['O2'], atom4= 'ID '+U_ID['C2'])) #dihedral angle: G.C6-G.N1-U.O2-U.C2
                    df.loc[((df['PDB_ID']==p1)&(df['bp_res']==x2)), "N1_O2_dih"]= N1_O2_dih
            
                elif n==2: #shifted G•U wobble 
                    #b1= N1-O4
                    #b1_dis= (cmd.get_distance(atom1= 'ID ' +G_ID['N1'], atom2= 'ID '+U_ID['O4'])) #distance between G.N1 and U.O4
                    N1_O4_dih= (cmd.get_dihedral(atom1= 'ID ' +G_ID['C2'], atom2= 'ID '+G_ID['N1'], atom3= 'ID '+U_ID['O4'], atom4= 'ID '+U_ID['C4'])) #dihedral angle: G.C2-G.N1-U.O4-U.C4
                    df.loc[((df['PDB_ID']==p1)&(df['bp_res']==x2)), "N1_O4_dih"]= N1_O4_dih
                    #b2= N2-N3
                    #b2_dis= (cmd.get_distance(atom1= 'ID ' +G_ID['N2'], atom2= 'ID '+U_ID['N3'])) #distance between G.N2 and U.N3
                    N2_N3_dih= (cmd.get_dihedral(atom1= 'ID ' +G_ID['C2'], atom2= 'ID '+G_ID['N2'], atom3= 'ID '+U_ID['N3'], atom4= 'ID '+U_ID['C4'])) #dihedral angle: G.C2-G.N1-U.N3-U.C4
                    df.loc[((df['PDB_ID']==p1)&(df['bp_res']==x2)), "N2_N3_dih"]= N2_N3_dih
            
                #print (b0_dis)
                #print (b0_dih)
                #print (b1_dis)
                #print (b1_dih)
                #print (b2_dis)
                #print (b2_dih)
                print (a2)
                print (b2)

                #extracting chain lengths for sequences in the structures
                #because in our study, both G and U are from the same chain
                #we are calculating the chain length only once
                stored.aseq={}
                cmd.iterate('chain '+a2['chain']+ ' and segi '+ a2['segi'], 'stored.aseq[resi]= resn')
                whole_aseq= ''
                for s1 in stored.aseq:
                    whole_aseq += stored.aseq[s1]

            
                #to save time we can skip the steps above to calculate temperature factors for the whole chain
                #in that case we can just run the following steps (mute the steps above to calculate the temperature factors for the whole chain) to claculate the temperature factors only for the G and U from the register of interest
                #for now we will just mute the following lines to calculate temperature factors calculation only for the the G and U from register of interest 
                stored.Gbs={}
                stored.Ubs={}
            
                for ids in G_ID:
                    if ids== "C1'":
                        pass
                    else:
                        cmd.iterate('chain '+ G_info['chain'] + ' and segi '+ G_info['segi'] + ' and resi '+ 'G' + ' and resn '+ G_info['resnum'] + 'and ID ' +G_ID[ids], 'stored.Gbs[ID]=b')
            
                for ids in U_ID:
                    if ids== "C1'":
                        pass
                    else:
                        cmd.iterate('chain '+ U_info['chain'] + ' and segi '+ U_info['segi'] + ' and resi '+ 'U' + ' and resn '+ U_info['resnum'] + 'and ID ' +U_ID[ids], 'stored.Ubs[ID]=b')

                G_b= sum([stored.Gbs[z3] for z3 in stored.Gbs])/len(stored.Gbs)
                U_b= sum([stored.Ubs[z3] for z3 in stored.Ubs])/len(stored.Ubs)
                b_comb= (G_b+U_b)/2

                #populating chain length for sequences in the structures
                df.loc[((df['PDB_ID']==p1)&(df['bp_res']==x2)), "Chain_length_structure"] = len(whole_aseq)

                #populating temperature factors for G, U, average of G and U, and whole molecule (structure)
                df.loc[((df['PDB_ID']==p1)&(df['bp_res']==x2)), "b_resG"] =  G_b
                df.loc[((df['PDB_ID']==p1)&(df['bp_res']==x2)), "b_resU"] =  U_b
                df.loc[((df['PDB_ID']==p1)&(df['bp_res']==x2)), "b_comb"] =  b_comb
                
                #notes on June 13, 2024
                print (distances)
 
            else:
                pass
   
        cmd.delete('all')
        q= p1.lower()+'.cif'
        os.remove(q) #removing the downloaded cif files
        
    return (df)

##optparser = OptionParser()
##(options, args) = optparser.parse_args()
##csvfile = args[0]



optparser = OptionParser()
(options, args) = optparser.parse_args()
csvfile = args[0]
D1= pd.read_csv(csvfile, keep_default_na=False, dtype={'seg_ID': str, 'res_ID_res1': str, 'seg_ID': str, 'res_ID_res2': str})

R= int(args[1])

D2= fix_PDB_ID(D1)

D3= get_info_structure(D2, R)

#store the finalized data as csv
#store the finalized data as csv
if R==1:
    #for standard wobble
    D3.to_csv('../results/all_standard_wobble_structure_data.csv')
elif R==2:
    #for shifted wobble
    D3.to_csv('../results/all_shifted_wobble_structure_data.csv')
#importing required libraries
import os
import pandas as pd
import numpy as np
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

#this function will clean the organism and RNA type information 
def clean_data(df1):
    #remove the examples with stapled RNAs
    df1= df1[df1['Molecule'].str.contains("stapled") == False]
    df1.index= np.arange(0, len(df1))
    print (df1.shape)
    
    #rename column --> 'Source_organism' to 'Source_organism_all' 
    df1.rename(columns = {'Source_Organism':'Source_Organism_all'}, inplace = True)

    
    #adding new columns--> Source_Organism_chain
    if 'Source_Organism_chain' in df1.columns:
        pass
    else:
        df1.insert (6, 'Source_Organism_chain', '')

    
    #fill out null data points with 'N/A'
    df1.fillna('N/A', inplace=True)
    
    #fill out 'unidentified' with 'N/A'
    df1['Source_Organism_all']= df1['Source_Organism_all'].replace('unidentified', 'N/A')
    
    df1['Expressed_Organism']= df1['Expressed_Organism'].replace('unidentified', 'N/A')

    
    #dictionary where key is old version of source organism and value is cleaned version of source organism
    org_cor={'Bacillus subtilis':'Bacillus subtilis', 'coli':'Escherichia coli',\
             'Thermus thermophilus':'Thermus thermophilus', 'Haloarcula marismortui': 'Haloarcula marismortui',\
             'Saccharomyces cerevisiae': 'Saccharomyces cerevisiae', 'Acinetobacter baumannii': 'Acinetobacter baumannii',\
             'Alicyclobacillus acidoterrestris': 'Alicyclobacillus acidoterrestris', 'Aquifex aeolicus': 'Aquifex aeolicus',\
             'Archaeoglobus fulgidus': 'Archaeoglobus fulgidus', 'Bombyx mori': 'Bombyx mori', \
             'Caldanaerobacter subterraneus subsp. tengcongensis': 'Caldanaerobacter subterraneus subsp. tengcongensis',\
             'Candida albicans':'Candida albicans', 'Deinococcus radiodurans':'Deinococcus radiodurans',  \
             'Desulfovibrio vulgaris':'Desulfovibrio vulgaris', 'Enterococcus faecalis':'Enterococcus faecalis',\
             'Fusobacterium nucleatum':'Fusobacterium nucleatum', 'Geobacillus stearothermophilus':'Geobacillus stearothermophilus',\
             'Halalkalibacterium halodurans':'Halalkalibacterium halodurans', 'Influenza A virus': 'Influenza A virus',\
             'Influenza B virus':'Influenza B virus', 'La Crosse': 'La Crosse virus',\
             'Lachnospiraceae bacterium':'Lachnospiraceae bacterium', 'Lactococcus lactis':'Lactococcus lactis',\
             'Listeria monocytogenes':'Listeria monocytogenes', 'Methanocaldococcus jannaschii':'Methanocaldococcus jannaschii',\
             'Methanocaldococcus jannaschii':'Methanocaldococcus jannaschii',\
             'Mycobacterium tuberculosis':'Mycobacterium tuberculosis',\
             'Neurospora crassa':'Neurospora crassa', 'Oceanobacillus iheyensis':'Oceanobacillus iheyensis',\
             'Pseudomonas aeruginosa':'Pseudomonas aeruginosa', 'Pyrococcus furiosus':'Pyrococcus furiosus',\
             'Salmonella enterica':'Salmonella enterica', 'Schizosaccharomyces pombe':'Schizosaccharomyces pombe',\
             'Staphylococcus aureus':'Staphylococcus aureus', 'Streptococcus sp.':'Streptococcus pyogenes',\
             'Mycobacterium smegmatis': 'Mycobacterium smegmatis',\
             'Leishmania donovani':'Leishmania donovani', 'Spinacia oleracea': 'Spinacia oleracea',\
             'Encephalitozoon cuniculi':'Encephalitozoon cuniculi',\
             'Homo sapiens':'Homo sapiens', 'Spinacia oleracea':'Spinacia oleracea',\
             'Cutibacterium acnes': 'Cutibacterium acnes', 'Spinacia oleracea':'Spinacia oleracea',\
             'HALOARCULA MARISMORTUI':'Haloarcula marismortui', 'Mus musculus':'Mus musculus',\
             'Drosophila melanogaster':'Drosophila melanogaster', 'Danio rerio': 'Danio rerio',
             'Oryctolagus cuniculus':'Oryctolagus cuniculus', 'Drosophila melanogaster':'Drosophila melanogaster',\
             'Xenopus laevis':'Xenopus laevis', 'Solanum lycopersicum':'Solanum lycopersicum',\
             'Cutibacterium acnes':'Cutibacterium acnes', 'Rattus norvegicus': 'Rattus norvegicus',\
             'Spraguea lophii':'Spraguea lophii', 'unidentified': 'unidentified',\
             'Thermococcus celer':'Thermococcus celer','THERMUS THERMOPHILUS': 'Thermus thermophilus',\
             'Streptococcus thermophilus':'Streptococcus thermophilus', 'Thermococcus kodakarensis':'Thermococcus kodakarensis',\
             'Thermococcus onnurineus':'Thermococcus onnurineus', 'Thermotoga maritima':'Thermotoga maritima',\
             'Tobacco mosaic virus':'Tobacco mosaic virus', 'Trypanosoma brucei':'Trypanosoma brucei',\
             'Vaccinia virus':'Vaccinia virus', 'Vibrio cholerae':'Vibrio cholerae', 'synthetic RNA': 'synthetic construct',\
            'Mycolicibacterium smegmatis':'Mycolicibacterium smegmatis'}
    
    #dictionary where key is old version of Molecule and value is cleaned version of Molecule
    mol_cor={'5S':'5S rRNA','5s': '5S rRNA', '16S':'16S rRNA', '16s':'16S rRNA', '23S': '23S rRNA', '23s':'23S rRNA',\
            '18S':'18S rRNA', '18s':'18S rRNA', '5.8S':'5.8S rRNA', '5.8s':'5.8S rRNA', '28S':'28S rRNA', '28s':'28S rRNA',\
            'tRNA':'tRNA', 'TRNA':'tRNA', 'messenger RNA':'mRNA', 'MRNA':'mRNA', 'transfer RNA':'tRNA', 'TRANSFER RNA':'tRNA',\
            '21S':'21S rRNA', '21s':'21S rRNA', '12s':'12S rRNA', '12S':'12S rRNA', '23S':'23S rRNA', '50S':'50S rRNA',\
            '16 ribosomal RNA':'16S rRNA', '2S ribosomal RNA':'2S rRNA', '23 ribosomal RNA':'23S rRNA', 'TRANSFER RNA':'tRNA'}
    
    #cleaning expressed organism
    for i, j in enumerate(df1.columns):
        if 'Expressed_Organism' in j:
            for k, l in enumerate(df1[j]):
                print (k)
                for n in org_cor:
                    if n in l:
                        print (l)
                        df1.loc[k, j]=org_cor[n]   
                        
    #cleaning source organism
    #source organism information which contain multiple organisms will not be cleaned
    for i, j in enumerate(df1.columns):
        if 'Source_Organism' in j:
            for k, l in enumerate(df1[j]):
                if ',' in l: #this indicates multiple organisms are present in source organism information
                    pass
                else:
                    for n in org_cor:
                        if n in l:
                            print (l)
                            df1.loc[k, j]=org_cor[n]
                  
    #cleaning Molecule (RNA type) information
    for i, j in enumerate(df1.columns):
        if 'Molecule' in j:
            for k, l in enumerate(df1[j]):
                if '25S' in l or '25s' in l:
                    df1.loc[k, j]= '#$%' #becasue '5S' is present as a key in mol_cor, '25S' cannot be corrected using as a key in mol_cor
                    #therefore, here we are giving it a unique name '#$%'
                    #in next loop other name will be corrected 
                    #finally '#$%' will be changed to '25S rRNA'
                else:
                    pass
    
    for i, j in enumerate(df1.columns):
        if 'Molecule' in j:
            for k, l in enumerate(df1[j]):
                for n in mol_cor:
                    if n in l:
                        df1.loc[k, j]=mol_cor[n]
    
    for i, j in enumerate(df1.columns):
        if 'Molecule' in j:
            for k, l in enumerate(df1[j]):
                if l=='#$%':
                    df1.loc[k, j] = '25S rRNA'
     
    #populating Source_Organism_chain
    for i, j in enumerate(df1['Expressed_Organism']):
        if j in df1['Source_Organism_all'][i]:
            df1.loc[i, 'Source_Organism_chain']= j
        elif j == 'N/A':
            df1.loc[i, 'Source_Organism_chain']= j
        else:
            df1.loc[i, 'Source_Organism_chain']= df1['Source_Organism_all'][i]
            
    df1['Source_Organism_chain'] = df1['Source_Organism_chain'].str.replace(',', '')
                             
                            
    return df1

#this function will filter the examples by applying different filtration criteria based on the registers of wobble
def filter_wobbles(df1, n):
    #filter out by distance, dihedral angle, and temperature factor cut-off
    #applying distance cut-off (for b1 and b2 -->less than of equal 3.4 A)
    #applying dihedral angle cut-off (for b1 and b2 --> less than or equal to positive or negative 50 degrees)
    #applying b0 distance cut-off --> distance between O6 of G and O4 of U has to be more than 3.4 A, therefore no interaction between these two atoms and we can confirm the formation of anionic reg 3
    #applying temperature factor cut-off --> excluding examples with temperature factors higher than 100 A2
    
    if n==1:
        b1_dis= 'O6_N3_dis'
        b1_dih= 'O6_N3_dih'
        b2_dis= 'N1_O2_dis'
        b2_dih= 'N1_O2_dih'
        
    
    elif n==2:
        b1_dis= 'N1_O4_dis'
        b1_dih= 'N1_O4_dih'
        b2_dis= 'N2_N3_dis'
        b2_dih= 'N2_N3_dih'
        

    #filter 1: applying b1, b2 distance cut-off
    df2= df1[(df1[b1_dis]<=3.4) & (df1[b2_dis]<=3.4)]
    df2.index = np.arange(0, len(df2))
    print (df2.shape)

    #filter 2: applying b1, b2 dihedral angle cut-off
    df3= df2[(abs(df2[b1_dih])<=50) & (abs(df2[b2_dih])<=50)]
    df3.index = np.arange(0, len(df3))
    print (df3.shape)

    #filter 3: 
    df4= df3[(df3['b_resG']<=100) & (df3['b_resU']<=100)]
    df4.index= np.arange(0, len(df4))
    print (df4.shape)
    
  
    #filter 4: applying b0 distance cut-off
    df5= df4[df4['O6_O4_dis'] > 3.4]
    df5.index= np.arange(0, len(df5))
    print (df5.shape)
    return df5

#read the data
optparser = OptionParser()
(options, args) = optparser.parse_args()
csvfile = args[0]
D1= pd.read_csv(csvfile)

R= int(args[1])

D2= fix_PDB_ID(D1)

#droping empty columns 
D2.dropna(inplace = True)
D2.index= np.arange(0, len(D2))

D3= filter_wobbles(D2, R)

D4= clean_data(D3)

if R==1:
    #for standard wobble
    D4.to_csv('../../results/all_red_good_standard_wobble.csv', index= False)
elif R==2:
    #for shifted wobble
    D4.to_csv('../../results/all_red_good_shifted_wobble.csv', index= False)


#importing required libraries
import os
import pandas as pd
import numpy as np
from optparse import OptionParser
from scipy.stats import percentileofscore

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


def cc_csv(dr, chain_ID):
    df1= pd.read_csv(dr)

    #converting multiple spaces to one space
    df1['Residue'] = df1['Residue'].str.replace(r'\s+', ' ', regex=True)
    df1['Residue'] = df1['Residue'].str.lstrip(' ') #removing one space from at the beginning of all residue entry
    df1.columns = df1.columns.str.lstrip(' ') #removing one space from at the beginning of all column name  

    #only keeping the chain of interest
    df2 = df1[df1['Residue'].str.startswith(chain_ID)]
    df2.index= np.arange(0, len(df2))

    #splitting the Residue (chain ID, residue index, residue ID, exactly in this order) column into two column: 'chain_ID', and 'residue_info'
    #'residue_info' contain the residue index and residue ID seperated by one space
    #this step is required, because in some rows there is no space between the chain_ID and residue index
    #however, there will always be space between residue index and residue ID
    df2[['chain_ID', 'residue_info']] = df2['Residue'].str.extract(f'^({chain_ID})(.*)', expand=True)

    #now splitting residue_info column into two columns: residue_index and residue_ID
    df2['residue_info'] = df2['residue_info'].str.lstrip(' ') 
    df2[['residue_index', 'residue_ID']] = df2['residue_info'].str.split(' ', expand=True)

    #updating the name for map-model correlation coefficient from ' CC' to residue_map_model_CC
    df2 = df2.rename(columns={'CC': 'residue_map_model_CC'})

    #keeping only the columns of interest
    df2 = df2.drop(columns= ['Residue', 'B_iso', 'Occupancy', '2Fo-Fc', 'Fmodel', 'residue_info'])
    df2 = df2[['chain_ID', 'residue_index', 'residue_ID','residue_map_model_CC' ]] #reordering columns 
    
    #rounding the CC values to 2 decimal places
    df2['residue_map_model_CC'] = df2['residue_map_model_CC'].astype(float)
    df2['residue_map_model_CC'] = df2['residue_map_model_CC'].round(2)

    #removing non-RNA residues, or any modified residues
    RNA_res= ['A', 'C', 'G', 'U']
    df3 = df2[df2['residue_ID'].isin(RNA_res)]
    df3.index= np.arange(0, len(df3))

    #removing any residues index with icode
    df4 = df3[~df3['residue_index'].str.contains(r'[a-zA-Z]', na=False)]

    #converting all residue_index into integer
    df4['residue_index'] = df4['residue_index'].astype(int)

    return df4


def cc_txt(dr, chain_ID):

    df1= pd.read_csv(dr)
    
    #converting multiple spaces to one space
    df1['CC per residue'] = df1['CC per residue'].str.replace(r'\s+', ' ', regex=True)
    
    #identify the location of rows from where 'CC per residue' information started
    ind_r = df1[df1['CC per residue'].str.contains('CC per chain')].index


    #################################
    #work here!!!!!!!!!!!!!!!!!!!!
    #removing the data for the 'CC per residue'
    df2 = df1.loc[:ind_r[0]-1]
    df2.index= np.arange(0, len(df2))
    
    #splitting the 'CC per residue' to five new columns
    #df3= pd.DataFrame(columns= ['A', 'chain_ID', 'residue_index', 'residue_ID', 'residue_map_model_CC'])
    
    df_check= df2['CC per residue'].str.split(' ', expand=True)
    if len(df_check.columns)==5:
        df3= pd.DataFrame(columns= ['A', 'chain_ID', 'residue_index', 'residue_ID', 'residue_map_model_CC'])
        df3[['A', 'chain_ID', 'residue_index', 'residue_ID', 'residue_map_model_CC']] = df2['CC per residue'].str.split(' ', expand=True)
        df3 = df3.drop(columns= ['A'])
    
    if len(df_check.columns)==4:
        df3= pd.DataFrame(columns= ['chain_ID', 'residue_index', 'residue_ID', 'residue_map_model_CC'])
        df3[['chain_ID', 'residue_index', 'residue_ID', 'residue_map_model_CC']] = df2['CC per residue'].str.split(' ', expand=True)
    
    #keeping only the chain of interest
    df4 = df3[df3['chain_ID']== chain_ID]
    df4.index= np.arange(0, len(df4))
    
    #rounding the CC values to 2 decimal places
    df4['residue_map_model_CC'] = df4['residue_map_model_CC'].astype(float)
    df4['residue_map_model_CC'] = df4['residue_map_model_CC'].round(2)

    #removing non-RNA residues, or any modified residues
    RNA_res= ['A', 'C', 'G', 'U']
    df5 = df4[df4['residue_ID'].isin(RNA_res)]
    df5.index= np.arange(0, len(df5))

    #removing any residues index with icode
    df6 = df5[~df5['residue_index'].str.contains(r'[a-zA-Z]', na=False)]

    #converting all residue_index into integer
    df6['residue_index'] = df6['residue_index'].astype(int)

    #for ind, item in enumerate(all_res_index):
    #    if ind== 0:
    #        good_res_index.append(item)
    #    else:
    #        if all_res_index[ind]- all_res_index[ind-1]< 1000:
    #            good_res_index.append(item)
    #        else:
    #            break

    #df2= df1[df1['residue_index'].isin(good_res_index)]
    #df2.index= np.arange(0, len(df2))

    #Removing MG
    #df1 = df1[df1['residue_ID'] != 'MG']
    #df1.index= np.arange(0, len(df1))

    #Removing CL
    #df1 = df1[df1['residue_ID'] != 'CL']
    #df1.index= np.arange(0, len(df1))

    #Removing NA
    #df1 = df1[df1['residue_ID'] != 'NA']
    #df1.index= np.arange(0, len(df1))

    #calculating z_score for G & U
    #df2['z_score'] = stats.zscore(df1['residue_map_model_CC'])
    #df2.z_score = df2.z_score.round(2)
    
    return df6


#read the data
optparser = OptionParser()
(options, args) = optparser.parse_args()
csvfile = args[0]
D1= pd.read_csv(csvfile)

cc_dir= args[1]

D2= fix_PDB_ID(D1)

sw_df1= D2[D2['representative']==1]
sw_df1.index= np.arange(0, len(sw_df1))

sw_df1['bp_res'] = sw_df1['bp_res'].str.replace('.', '_', regex=False)
sw_df1['bp_res'] = sw_df1['bp_res'].str.replace('-', '_', regex=False)

home= os.getcwd()

os.chdir(cc_dir)

stat_count=0
for i, j in enumerate(sw_df1['PDB_ID']):
    stat_count +=1
    print ('working here----> '+str(stat_count))
    print (j)
    t_fname= j+ ".txt"
    l_fname= j+ ".log"
    c_fname= j+ ".csv"
    if os.path.isfile(t_fname):
        print ('#################################')
        print (sw_df1['chain_ID'][i])
        df1= cc_txt(t_fname, sw_df1['chain_ID'][i])
        #df1= run_cryoEM(df)
        #print (df1)
    #elif os.path.isfile(l_fname):
    #    print ('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
    #    df1= run_xray(l_fname)
        #print (df1)
    
    elif os.path.isfile(c_fname):
        print ('#################################')
        df1= cc_csv(c_fname, sw_df1['chain_ID'][i])

    
    print (df1.shape)
    #filtering out the chain of interest
    #df2= df1[df1['chain_ID']== sw_df1['chain_ID'][i]]
    #df2.index= np.arange(0, len(df2))

    #residue index for G
    if sw_df1['res_ID_res1'][i]== 'G':
        res_index_G= sw_df1['res_index_res1'][i]
        res_index_U= sw_df1['res_index_res2'][i]
    else:
        res_index_G= sw_df1['res_index_res2'][i]
        res_index_U= sw_df1['res_index_res1'][i]

    cc_median = round(df1['residue_map_model_CC'].median(), 2)
    cc_mean = round(df1['residue_map_model_CC'].mean(),2)

    cc_G= df1[df1['residue_index']== res_index_G]['residue_map_model_CC'].to_list()[0]
    cc_U= df1[df1['residue_index']== res_index_U]['residue_map_model_CC'].to_list()[0]

    percentile_G = percentileofscore(df1['residue_map_model_CC'], cc_G)
    percentile_U = percentileofscore(df1['residue_map_model_CC'], cc_U)

    D2.loc[((D2['PDB_ID']== j)&(D2['chain_ID']== sw_df1['chain_ID'][i])), 'Median_map_model_cc'] = cc_median
    D2.loc[((D2['PDB_ID']== j)&(D2['chain_ID']== sw_df1['chain_ID'][i])), 'Mean_map_model_cc'] = cc_mean

    D2.loc[((D2['PDB_ID']== j)&(D2['chain_ID']== sw_df1['chain_ID'][i])), 'Raw_map_model_cc_for_the_G'] = cc_G
    D2.loc[((D2['PDB_ID']== j)&(D2['chain_ID']== sw_df1['chain_ID'][i])), 'Normalized_map_model_cc_for_the_G'] = round(cc_G/cc_mean, 2)
    D2.loc[((D2['PDB_ID']== j)&(D2['chain_ID']== sw_df1['chain_ID'][i])), 'Percentile_cc_for_the_G'] = round(percentile_G, 2)

    D2.loc[((D2['PDB_ID']== j)&(D2['chain_ID']== sw_df1['chain_ID'][i])), 'Raw_map_model_cc_for_the_U'] = cc_U
    D2.loc[((D2['PDB_ID']== j)&(D2['chain_ID']== sw_df1['chain_ID'][i])), 'Normalized_map_model_cc_for_the_U'] = round(cc_U/cc_mean, 2)
    D2.loc[((D2['PDB_ID']== j)&(D2['chain_ID']== sw_df1['chain_ID'][i])), 'Percentile_cc_for_the_U'] = round(percentile_U, 2)

os.chdir(home)

D2.to_csv('data/nr_shifted_wobble_map_model_cc.csv', index= False)
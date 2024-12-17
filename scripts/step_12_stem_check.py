#importing required libraries
import os
from subprocess import call
import pandas as pd
import json
import numpy as np
from optparse import OptionParser

def get_bps(jfile, chain_ID, res1_index, res2_index):
    #this function will extract the base pair of interest for five locations (two base pair next to the wobble, wobble base pair itself, and two base pair before the wobble)
    with open(jfile, 'r') as f:
        data= json.loads(f.read())
    bp_list= pd.json_normalize(data, record_path =['pairs']) #storing base pair information into a dataframe

    #inserting additional column for nt 1
    bp_list.insert(2, 'nt1_chain', bp_list['nt1'].astype(str).str.split('.', expand=True)[0])
    bp_list.insert(3, 'nt1_resi', bp_list['nt1'].astype(str).str.split('.', expand=True)[1].astype(str).str[1:])
    
    #inserting additional column for nt 2
    bp_list.insert(5, 'nt2_chain', bp_list['nt2'].astype(str).str.split('.', expand=True)[0])
    bp_list.insert(6, 'nt2_resi', bp_list['nt2'].astype(str).str.split('.', expand=True)[1].astype(str).str[1:])

    #applying chain ID filter
    #this function is only considering base pairs within the same chain
    bp_list_1= bp_list[(bp_list["nt1_chain"]==chain_ID)& (bp_list["nt2_chain"]==chain_ID)]
    bp_list_1.index = np.arange(0, len(bp_list_1))

    #applying base pair name filter
    bp_list_2= bp_list_1[(bp_list_1['name']=='WC')|(bp_list_1['nt1_resi']== str(res1_index))]
    bp_list_2.index = np.arange(0, len(bp_list_2))

    bp_list_3= bp_list_2[(bp_list_2['name']=='WC')|(bp_list_2['nt2_resi']== str(res2_index))]
    bp_list_3.index = np.arange(0, len(bp_list_3))
    row_index = bp_list_3.index[bp_list_3['nt1_resi'] == str(res1_index)].to_list()

    ind1= row_index[0]-2
    ind3= row_index[0]+2
    bp_list_4= bp_list_3.loc[ind1: ind3]
    bp_list_4.index = np.arange(0, len(bp_list_4))

    return bp_list_4


#read the csv file
optparser = OptionParser()
(options, args) = optparser.parse_args()
csvfile = args[0]
D1= pd.read_csv(csvfile, keep_default_na=False, dtype={'seg_ID': str, 'res_ID_res1': str, 'seg_ID': str, 'res_ID_res2': str})

#get the directory containing all json files 
json_dr = args[1]

R= int(args[2]) #R is indicating the register we are analyzing here, 1= standard wobble, 2= shifted wobble

#below is the home directory, this script is srtored in this directory
hm_dr= os.getcwd()


count =0
os.chdir(json_dr)
sec_str_loc= {'unstructured': [], 'stem': [], 'unpaired_above_GU': [], 'unpaired_above_UG': [], 'unpaired_below_GU': [], 'unpaired_below_UG':[], 'inside_loop':[], 'unidentified': []}

D1['bp_1_above']= ''
D1['bp_2_above']= ''
D1['bp_1_below']= ''
D1['bp_2_below']= ''

for i, j in enumerate(D1['bp_ID']):
    count +=1
    print (count)
    print (str(count)+'. working here ------------------------------------------> '+ str(j))
    p_ID= j.split('_')[0]
    #7M4Y_a_U1083-a_G1096
    res1= j.split('-')[0].split('_')[2]
    res1_index= int(res1[1:]) #index for residue 1
    res1_above= res1_index+1 #index for residue one step above of residue 1
    res1_above1= res1_index+2#index for residue two steps above of residue 1
    res1_below= res1_index-1 #index for residue one step below of residue 1
    res1_below1= res1_index-2#index for residue two steps below of residue 1

    res2= j.split('-')[1].split('_')[1]
    res2_index= int(res2[1:]) #index for residue 2
    res2_above= res2_index-1 #index for residue one step above of residue 2
    res2_above1= res2_index-2#index for residue two steps above of residue 2
    res2_below= res2_index+1 #index for residue one step below of residue 2
    res2_below1= res2_index+2#index for residue two steps below of residue 2
    
    chain_ID= j.split('_')[1] #chain ID is same for both res1 and res2
    
    jfile= p_ID+'-dssr.json' #name of the corresponding JSON file, this file contain all characterization output
    
    df_bps= get_bps(jfile, chain_ID, res1_index, res2_index)
    print (df_bps)
    
    df_above= df_bps[(df_bps['nt1_resi']== str(res1_above))&(df_bps['nt2_resi']== str(res2_above))] #bp one step above
    df_above1= df_bps[(df_bps['nt1_resi']== str(res1_above1))&(df_bps['nt2_resi']== str(res2_above1))] #bp two steps above

    df_below= df_bps[(df_bps['nt1_resi']== str(res1_below))&(df_bps['nt2_resi']== str(res2_below))] #bp one step below
    df_below1= df_bps[(df_bps['nt1_resi']== str(res1_below1))&(df_bps['nt2_resi']== str(res2_below1))] #bp two steps below

    if len(df_above)>0:
        bp_above= df_above['bp'].to_list()[0]
        D1.loc[(D1['bp_ID']==j), "bp_1_above"] =  bp_above #bp one step above

    if len(df_above1)>0:
        bp_above1= df_above1['bp'].to_list()[0]
        D1.loc[(D1['bp_ID']==j), "bp_2_above"] =  bp_above1 #bp two steps above

    if len(df_below)>0:
        bp_below= df_below['bp'].to_list()[0]
        D1.loc[(D1['bp_ID']==j), "bp_1_below"] =  bp_below #bp one step below

    if len(df_below1)>0:
        bp_below1= df_below1['bp'].to_list()[0]
        D1.loc[(D1['bp_ID']==j), "bp_2_below"] =  bp_below1 #bp two steps below

    

    

    if len(df_bps)== 5:
        #df_above= df_bps[(df_bps['nt1_resi']== str(res1_above))&(df_bps['nt2_resi']== str(res2_above))]
        #########df_above1

        #df_below= df_bps[(df_bps['nt1_resi']== str(res1_below))&(df_bps['nt2_resi']== str(res2_below))]
        #########df_below1
        flnk_ind= []
        flnk_ind.append(df_bps.loc[3]['nt1_resi'])
        flnk_ind.append(df_bps.loc[1]['nt1_resi'])
        flnk_ind.append(df_bps.loc[3]['nt2_resi'])
        flnk_ind.append(df_bps.loc[1]['nt2_resi'])
        flnk_ind_1= []

        for item in flnk_ind:
            if not str(item).isdigit():
                pass
            else:
                flnk_ind_1.append(item)
    
        if len(flnk_ind_1)==4:
            left_top= int(flnk_ind_1[0])
            left_bottom= int(flnk_ind_1[1])
            right_top= int(flnk_ind_1[2])
            right_bottom= int(flnk_ind_1[3])




        #left_top= int(df_bps.loc[3]['nt1_resi'])
        #left_bottom= int(df_bps.loc[1]['nt1_resi'])

        #right_top= int(df_bps.loc[3]['nt2_resi'])
        #right_bottom= int(df_bps.loc[1]['nt2_resi'])
    
            if len(df_below)>0 and len(df_above)>0:
                #these are examples where wobbles are within a stem
                sec_str_loc['stem'].append(j)

                #bp_above= df_above['bp'].to_list()[0]
                #print (bp_above)
                #bp_below= df_below['bp'].to_list()[0]
                #print (bp_below)
                #print (df_bps)

                #D1.loc[(D1['bp_ID']==j), "bp_above"] =  bp_above
                #########D1.loc[(D1['bp_ID']==j), "bp_above1"] =  bp_above1
                #D1.loc[(D1['bp_ID']==j), "bp_below"] =  bp_below
                #########D1.loc[(D1['bp_ID']==j), "bp_below1"] =  bp_below1
    
            if len(df_above)==0 and len(df_below)==0:
        
                if abs(left_top-left_bottom)<11 and abs(right_top-right_bottom)<11:
                    #these are examples where wobbles are within a internal loops or bulges
                    sec_str_loc['inside_loop'].append(j)
            
                elif abs(right_bottom-left_bottom)<11:
                    #these are examples where wobbles are within a hairpin loops
                    sec_str_loc['inside_loop'].append(j)
            
                else:
                    #these are examples where wobbles are within unstructured motifs  
                    sec_str_loc['unstructured'].append(j)

            if len(df_above)>0 and len(df_below)==0:
                #bp_above= df_above['bp'].to_list()[0]
                #D1.loc[(D1['bp_ID']==j), "bp_above"] =  bp_above
                #print (bp_above)
                #print (df_bps)
                #these are examples of terminal GU/UG where WCF base pairs are below the GU/UG
                if D1['res_ID_res1'][i] =='G':
                    #GU examples
                    sec_str_loc['unpaired_below_GU'].append(j)
                else:
                    #UG examples 
                    sec_str_loc['unpaired_below_UG'].append(j)
            
            
            if len(df_below)>0 and len(df_above)==0:
                #bp_below= df_below['bp'].to_list()[0]
                #D1.loc[(D1['bp_ID']==j), "bp_below"] =  bp_below
                #print (bp_below)
                #print (df_bps)
                #these are examples of terminal GU/UG where WCF pairs are above rthe GU/UG
                if (left_top in range(res1_index, res2_index)) and (right_top in range(res1_index, res2_index)):
                    pass
                if D1['res_ID_res1'][i] =='G':
                    #GU examples 
                    sec_str_loc['unpaired_above_GU'].append(j)
                else:
                    #UG examples
                    sec_str_loc['unpaired_above_UG'].append(j)
        else: 
            sec_str_loc['unidentified'].append(j)
    else:
        sec_str_loc['unidentified'].append(j)
os.chdir(hm_dr)   

def get_key(value):
    for key, values in sec_str_loc.items():
        if value in values:
            return key
D1['location_in_secondary_structure'] = D1['bp_ID'].apply(get_key)

if R==1:
    D1.to_csv('nr_standard_wobble_2D_loc.csv', index= False)
elif R==2:
    D1.to_csv('nr_shifted_wobble_2D_loc.csv', index= False)
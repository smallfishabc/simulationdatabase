import os
import numpy as np
import pandas as pd
import datacleaning
import database
import Generatedatabase
# Here, we will define several database operation that are frequently used in my previous scripts.

# Select concentration we want
def select_MTFE_data(data, MTFE):
    MTFEdata = data[data['MTFE'] == MTFE]
    return MTFEdata

# Select the data type we want (Rg,Re,...) Default concentration will be 0 (water)
def select_datatype_data(data, datatype, MTFE=0):
    typedata = data[data['datatype'] == datatype].reset_index()
    return typedata

# Select the protein we want
def select_protein_data(data, protein):
    protein_data = data[data['Protein'] == protein].reset_index()
    return protein_data

# Make several selection at the same time (Will modify this according to my scripts)
def multiple_selection_function(data,protein,datatype,MTFE):
    df1=select_protein_data(data,protein)
    df2=select_datatype_data(df1,datatype)
    df3=select_MTFE_data(df2,MTFE)
    return df3

# Calculate the Re of ideal polymer
def GS_length_calculation(sequence):
    length=len(sequence)
    return(0.55359657*(length)**0.47499067)

# Calculate the Re Ratio between target protein and ideal polymer
def chi_value_calculation(data,protein,MTFE,sequence):
    ee=multiple_selection_function(data,protein,'ee',MTFE)['Rs'].tolist()[0]
    l_GS=GS_length_calculation(sequence)
    chi=ee/l_GS-1
    return(chi)

# Calculate the solution sensitivity of proteins in attractive solutions
def attractive_sensitivity_calculation(data,protein,sequence,chi_0,datatype='ee',MTFE=3):
    chi_p3=chi_value_calculation(data,protein,MTFE,sequence)
    sensitivity= chi_p3-chi_0
    return sensitivity

# Calculate the solution sensitivity of proteins in repulsive solutions
def repulsive_sensitivity_calculation(data,protein,sequence,chi_0,datatype='ee',MTFE=-3):
    chi_m3=chi_value_calculation(data,protein,MTFE,sequence)
    sensitivity= chi_m3-chi_0
    return sensitivity

# A plot function designed for finding correlation between different datatype.
def database_plot_pre(datatype,MTFE=0):
    df_entry = pd.read_csv('database_entry.csv')
    df = pd.read_csv('database_full_value_test1.csv')
    full_df= pd.DataFrame()
    for i in range(len(df_entry)):
        protein_name = df_entry.loc[i, 'Protein']
        sequence = df_entry.loc[i, 'Sequence']
        valueraw1 = multiple_selection_function(df,protein_name,datatype,MTFE)['att1'].tolist()[0]
        valueraw2 = multiple_selection_function(df, protein_name, datatype, MTFE)['att2'].tolist()[0]
        valueraw3 = multiple_selection_function(df, protein_name, datatype, MTFE)['rep1'].tolist()[0]
        valueraw4 = multiple_selection_function(df, protein_name, datatype, MTFE)['rep2'].tolist()[0]
        # if valueraw==0:
        #     continue
        #value=(valueraw1+2*valueraw2-valueraw3-2*valueraw4)/len(sequence)
        #value = (valueraw1 + 2 * valueraw2) / len(sequence)
        #value = (valueraw1 + valueraw2) / len(sequence)
        value = (valueraw1 +  valueraw2 + valueraw3 +  valueraw4) / len(sequence)
        chi_0=chi_value_calculation(df,protein_name, 0,sequence)
        att_sensitivity=attractive_sensitivity_calculation(df,protein_name,sequence,chi_0)
        rep_sensitivity=repulsive_sensitivity_calculation(df, protein_name, sequence, chi_0)
        dict_pd={'Protein': protein_name,datatype:value,'chi_0':chi_0,'att_sensitivity':att_sensitivity,'rep_sensitivity':rep_sensitivity}
        new_data=pd.DataFrame(dict_pd,index=[0])
        full_df = pd.concat([full_df, new_data], ignore_index=True)
    return full_df
# Designed for calculating several value in on function
def database_plot_pre_multiple(datatype,MTFE=0,index=8):
    #os.chdir('/media/lemoncatboy/WD_BLACK/DATA_F/puma_scramble_new/puma_scrammble_sum')
    df_entry = pd.read_csv('database_entry.csv')
    df = pd.read_csv('database_full_value_0428_entropy.csv')
    #df = pd.read_csv('database_full_value_1021_iupred.csv')
    #df = pd.read_csv('database_full_value_1021_urgent.csv')
    #df = pd.read_csv('database_full_value_1020.csv')
    #df = pd.read_csv('database_full_value_1104_PDBsumreal.csv')
    full_df= pd.DataFrame()
    for i in range(len(df_entry)):
        protein_name = df_entry.loc[i, 'Protein']
        sequence = df_entry.loc[i, 'Sequence']
        seq_length = len(sequence)
        seq_length_rt=np.sqrt(seq_length)
        #iupred = df_entry.loc[i,'iupred']
        # Add the df_protein for selecting sequence feature
        df_protein = df[df['Protein'] == protein_name]
        # Covert list from csv to python list
        valueraw1_list, valueraw2_list, valueraw3_list, valueraw4_list=interaction_list_generate(df, protein_name, datatype, MTFE)
        valueraw1_list_p3, valueraw2_list_p3, valueraw3_list_p3, valueraw4_list_p3 = interaction_list_generate(df, protein_name,
                                                                                                   datatype, MTFE=3)
        valueraw1_list_m3, valueraw2_list_m3, valueraw3_list_m3, valueraw4_list_m3 = interaction_list_generate(df, protein_name,
                                                                                                   datatype, MTFE=-3)
        # Record the datatype of each value using typelist
        typelist=['none_none','none_ratio','none_probability','dis_none','dis_ratio','dis_probability','rtdis_none','rtdis_ratio','rtdis_probability']
        # Select data from the list
        type=typelist[index]

        interaction_data=interaction_pre_list(index, seq_length,typelist, valueraw1_list, valueraw2_list, valueraw3_list, valueraw4_list)
        interaction_data_p3 = interaction_pre_list(index, seq_length,typelist, valueraw1_list_p3, valueraw2_list_p3, valueraw3_list_p3,
                                                valueraw4_list_p3,MTFE=3)
        interaction_data_m3 = interaction_pre_list(index, seq_length,typelist, valueraw1_list_m3, valueraw2_list_m3, valueraw3_list_m3,
                                                valueraw4_list_m3,MTFE=-3)
        chi_0=chi_value_calculation(df,protein_name, 0,sequence)
        att_sensitivity=attractive_sensitivity_calculation(df,protein_name,sequence,chi_0,MTFE=3)
        rep_sensitivity=repulsive_sensitivity_calculation(df, protein_name, sequence, chi_0,MTFE=-3)
        heli=multiple_selection_function(df, protein_name, 'helix', 0)['Rs'].tolist()[0]
        HB = multiple_selection_function(df, protein_name, 'HB', 0)['Rs'].tolist()[0]
        beta = multiple_selection_function(df, protein_name, 'beta', 0)['Rs'].tolist()[0]
        # Calculate Kappa value
        delta = select_datatype_data(df_protein, 'feature', MTFE=0)['delta']
        #dict_pd={'Protein': protein_name,type+'att':value,type+'rep':value2,'chi_0':chi_0,'att_sensitivity':att_sensitivity,
        #       'rep_sensitivity':rep_sensitivity,'Helicity':heli,'H-bonds':HB,'delta':delta}
        dict_pd={'Protein': protein_name,'length':seq_length,'chi_0':chi_0,'att_sensitivity':att_sensitivity,
               'rep_sensitivity':rep_sensitivity,'Helicity':heli,'H-bonds':HB,'delta':delta,'Beta_sheet':beta}
        new_data = pd.DataFrame(dict_pd,index=[0])
        new_data = pd.concat([new_data,interaction_data,interaction_data_m3,interaction_data_p3],axis=1)
        full_df = pd.concat([full_df, new_data], ignore_index=True)
        full_df.to_csv('interaction_strength_fitting_1109.csv',index=False)
    return full_df
def plot_pre_list(datatype):
    list_data=[]
    for i in range(9):
        j=database_plot_pre_multiple(datatype,MTFE=0,index=i)
        list_data.append(j)
    return list_data
# To simplify the interaction calculation process, we created two separate functions
def interaction_list_generate(df, protein_name, datatype, MTFE):
    valueraw1_raw = multiple_selection_function(df, protein_name, datatype, MTFE)['att1'].tolist()[0]
    valueraw1_list = datacleaning.covert_string_tolist(valueraw1_raw)
    valueraw2_raw = multiple_selection_function(df, protein_name, datatype, MTFE)['att2'].tolist()[0]
    valueraw2_list = datacleaning.covert_string_tolist(valueraw2_raw)
    valueraw3_raw = multiple_selection_function(df, protein_name, datatype, MTFE)['rep1'].tolist()[0]
    valueraw3_list = datacleaning.covert_string_tolist(valueraw3_raw)
    valueraw3_list = datacleaning.absolute_value(valueraw3_list)
    valueraw4_raw = multiple_selection_function(df, protein_name, datatype, MTFE)['rep2'].tolist()[0]
    valueraw4_list = datacleaning.covert_string_tolist(valueraw4_raw)
    valueraw4_list = datacleaning.absolute_value(valueraw4_list)
    return (valueraw1_list,valueraw2_list,valueraw3_list,valueraw4_list)
def interaction_pre_list(index,seq_length,typelist,valueraw1_list,valueraw2_list,valueraw3_list,valueraw4_list,MTFE=0):
    full_data=pd.DataFrame({'a':1},index=[0])
    for i in range(index):
        valueraw1=valueraw1_list[i]
        valueraw2=valueraw2_list[i]
        valueraw3=valueraw3_list[i]
        valueraw4=valueraw4_list[i]
        #rawatt=(valueraw1+valueraw2)/seq_length
        rawatt = (valueraw1 + valueraw2)
        #rawrep=(valueraw3+valueraw4)/seq_length
        rawrep = (valueraw3 + valueraw4)
        column_name=typelist[i]
        column_name_att=column_name+'att'+str(MTFE)
        column_name_rep=column_name+'rep'+str(MTFE)
        full_data[column_name_att]=rawatt
        full_data[column_name_rep]=rawrep
    full_data.drop(columns=['a'],inplace=True)
    return full_data
def database_plot_pre_multiple_none_interaction(datatype,MTFE=0,index=8):
    #os.chdir('/media/lemoncatboy/WD_BLACK/DATA_F/puma_scramble_new/puma_scrammble_sum')
    df_entry = pd.read_csv('database_entry.csv')
    df = pd.read_csv('database_full_value_0428_entropy.csv')
    #df = pd.read_csv('database_full_value_1021_iupred.csv')
    #df = pd.read_csv('database_full_value_1021_urgent.csv')
    #df = pd.read_csv('database_full_value_1020.csv')
    #df = pd.read_csv('database_full_value_1104_PDBsumreal.csv')
    full_df= pd.DataFrame()
    for i in range(len(df_entry)):
        protein_name = df_entry.loc[i, 'Protein']
        sequence = df_entry.loc[i, 'Sequence']
        seq_length = len(sequence)
        seq_length_rt=np.sqrt(seq_length)
        #iupred = df_entry.loc[i,'iupred']
        # Add the df_protein for selecting sequence feature
        df_protein = df[df['Protein'] == protein_name]
        chi_0=chi_value_calculation(df,protein_name, 0,sequence)
        att_sensitivity=attractive_sensitivity_calculation(df,protein_name,sequence,chi_0,MTFE=3)
        rep_sensitivity=repulsive_sensitivity_calculation(df, protein_name, sequence, chi_0,MTFE=-3)
        heli=multiple_selection_function(df, protein_name, 'helix', 0)['Rs'].tolist()[0]
        HB = multiple_selection_function(df, protein_name, 'HB', 0)['Rs'].tolist()[0]
        beta = multiple_selection_function(df, protein_name, 'beta', 0)['Rs'].tolist()[0]
        # Calculate Kappa value
        #dict_pd={'Protein': protein_name,type+'att':value,type+'rep':value2,'chi_0':chi_0,'att_sensitivity':att_sensitivity,
        #       'rep_sensitivity':rep_sensitivity,'Helicity':heli,'H-bonds':HB,'delta':delta}
        dict_pd={'Protein': protein_name,'length':seq_length,'chi_0':chi_0,'att_sensitivity':att_sensitivity,
               'rep_sensitivity':rep_sensitivity,'Helicity':heli,'H-bonds':HB,'Beta_sheet':beta}
        new_data = pd.DataFrame(dict_pd,index=[0])
        full_df = pd.concat([full_df, new_data], ignore_index=True)
        full_df.to_csv('interaction_strength_fitting_1109.csv',index=False)
    return full_df
if __name__ == '__main__':
    directory = 'F:\DATA_F\GSlinker_entropic_force'
    os.chdir(directory)
    database_plot_pre_multiple_none_interaction('interaction')
import pandas as pd

import database
import Generatedatabase
# Here, we will define several database operation that are frequently used in my previous scripts.

def select_MTFE_data(data, MTFE):
    MTFEdata = data[data['MTFE'] == MTFE]
#    MTFEdata.to_csv('MTFE' + str(MTFE) + 'data.csv', index=False)
    return MTFEdata


def select_datatype_data(data, datatype, MTFE=0):
    typedata = data[data['datatype'] == datatype].reset_index()
#    typedata.to_csv('type' + str(datatype) + 'MTFE' + str(MTFE) + 'data.csv', index=False)
    return typedata


def select_protein_data(data, protein):
    protein_data = data[data['Protein'] == protein].reset_index()
    return protein_data

def multiple_selection_function(data,protein,datatype,MTFE):
    df1=select_protein_data(data,protein)
    df2=select_datatype_data(df1,datatype)
    df3=select_MTFE_data(df2,MTFE)
    return df3

def GS_length_calculation(sequence):
    length=len(sequence)
    return(0.55359657*(length)**0.47499067)

def chi_value_calculation(data,protein,MTFE,sequence):
    ee=multiple_selection_function(data,protein,'ee',MTFE)['Rs'].tolist()[0]
    print(ee)
    l_GS=GS_length_calculation(sequence)
    chi=ee/l_GS-1
    return(chi)

def attractive_sensitivity_calculation(data,protein,sequence,chi_0,datatype='ee'):
    chi_p3=chi_value_calculation(data,protein,3,sequence)
    sensitivity= chi_p3-chi_0
    return sensitivity

def repulsive_sensitivity_calculation(data,protein,sequence,chi_0,datatype='ee'):
    chi_m3=chi_value_calculation(data,protein,-3,sequence)
    sensitivity= chi_m3-chi_0
    return sensitivity

def database_plot_pre(datatype,MTFE=0):
    df_entry = pd.read_csv('database_entry.csv')
    df = pd.read_csv('database_full_value_test2.csv')
    full_df= pd.DataFrame()
    for i in range(len(df_entry)):
        protein_name = df_entry.loc[i, 'Protein']
        sequence = df_entry.loc[i, 'Sequence']
        valueraw1=multiple_selection_function(df,protein_name,datatype,MTFE)['att1'].tolist()[0]
        valueraw2 = multiple_selection_function(df, protein_name, datatype, MTFE)['att2'].tolist()[0]
        valueraw3 = multiple_selection_function(df, protein_name, datatype, MTFE)['rep1'].tolist()[0]
        valueraw4 = multiple_selection_function(df, protein_name, datatype, MTFE)['rep2'].tolist()[0]
        # if valueraw==0:
        #     continue
        #value=(valueraw1+2*valueraw2-valueraw3-2*valueraw4)/len(sequence)
        #value = (valueraw1 + 2 * valueraw2) / len(sequence)
        #value = (valueraw1 + valueraw2) / len(sequence)
        value = (valueraw1 +  valueraw2 + valueraw3 +  valueraw4) / len(sequence)
        print(value)
        chi_0=chi_value_calculation(df,protein_name, 0,sequence)
        att_sensitivity=attractive_sensitivity_calculation(df,protein_name,sequence,chi_0)
        rep_sensitivity=repulsive_sensitivity_calculation(df, protein_name, sequence, chi_0)
        dict_pd={'Protein': protein_name,datatype:value,'chi_0':chi_0,'att_sensitivity':att_sensitivity,'rep_sensitivity':rep_sensitivity}
        new_data=pd.DataFrame(dict_pd,index=[0])
        full_df = pd.concat([full_df, new_data], ignore_index=True)
    return full_df
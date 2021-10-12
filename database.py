import pandas as pd
import os
import pysnooper

def load_data(df):
    df_rg = add_data('BBrawRG_5en.csv', df, 'rg')
    df_ee = add_data('BBEERG_5en.csv', df, 'ee')
    df_HB = add_data_HB('BB_HB_5en.csv', df, 'HB')
    df_helix = add_data('BB_Heli_5en.csv', df, 'helix')
    df_interaction = add_data('BBcontact_new_test2.csv', df, 'interaction')
    df_feature = add_data_sequence_feature(df)
    #df_interaction = add_data('BB_contact_lines.csv',df,'interaction')
    #   df.contact = add_data('BB_Heli_5en.csv',df,'contact')
    fulldata = pd.concat([df_ee, df_rg, df_HB, df_helix,df_interaction,df_feature])
    #fulldata.to_csv('fulldata.csv', index=False)
    return fulldata
    # bufferdata = fulldata[fulldata['MTFE'] == 0.0]
    # bufferdata.to_csv('bufferdata.csv',index=False)
    # eedata = bufferdata[bufferdata['datatype'] == 'ee'].reset_index()
    # eedata.insert(eedata.shape[1], 'ratio', 'Nan')
    # rgdata = bufferdata[bufferdata['datatype'] == 'rg'].reset_index()
    # rgdata.to_csv('rgdata.csv',index=False)
    # for i in range(len(eedata)):
    #     protein_name = eedata.loc[i, 'Protein']
    #     print(protein_name)
    #     ee = eedata.loc[i, 'Rs']
    #     rg = rgdata[rgdata['Protein'] == protein_name]['Rs']
    #     print(ee)
    #     print(rg.iloc[0])
    #     ratio = ee / rg.iloc[0]
    #     eedata.loc[i, 'ratio'] = ratio
    # eedata.to_csv('eedata.csv',index=False)


def add_data(filename, entrydf, data_type):
    full_df = pd.DataFrame()
    for i in range(len(entrydf)):
        protein_name = entrydf.loc[i, 'Protein']
        protein_directory = entrydf.loc[i, 'Directory']
        new_data = read_data(protein_directory, filename)
        new_data.insert(0, 'datatype', data_type)
        new_data.insert(0, 'Protein', protein_name)
        full_df = pd.concat([full_df, new_data], ignore_index=True)
    return full_df
# Because the HB file has different format, we need to rename the HB when attach to the database

def add_data_HB(filename, entrydf, data_type):
    full_df = pd.DataFrame()
    for i in range(len(entrydf)):
        protein_name = entrydf.loc[i, 'Protein']
        protein_directory = entrydf.loc[i, 'Directory']
        new_data = read_data(protein_directory, filename)
        new_data.insert(0, 'datatype', data_type)
        new_data.insert(0, 'Protein', protein_name)
        full_df = pd.concat([full_df, new_data], ignore_index=True)
    full_df=full_df.rename(columns={'Hbond': 'Rs', 'st': 'Sd'}, errors="raise")
    return full_df

def add_data_sequence_feature(entrydf, data_type='feature'):
    full_df = pd.DataFrame()
    for i in range(len(entrydf)):
        protein_name = entrydf.loc[i, 'Protein']
        protein_directory = entrydf.loc[i, 'Directory']
        sequence= entrydf.loc[i,'Sequence']
        length=len(sequence)
        print(sequence,length)
        #Will calculate specific sequence feature
        #cider_calculation()
        #Feng_calculation()
        new_data=pd.DataFrame({'Protein':protein_name,'datatype':data_type,'Sequence':sequence,'length':length},index=[0])

        #new_data.insert(0, 'datatype', data_type)
        #new_data.insert(0, 'Protein', protein_name)
        full_df = pd.concat([full_df, new_data], ignore_index=True)
    return full_df


def read_data(directory, filename):
    filepath = os.path.join(directory, filename)
    data = pd.read_csv(filepath)
    return (data)



import pandas as pd
import os


def load_data(df):
    df_rg = add_data('BBrawRG_5en.csv', df, 'rg')
    df_ee = add_data('BBEERG_5en.csv', df, 'ee')
    df_HB = add_data('BB_HB_5en.csv', df, 'contact')
    df_helix = add_data('BB_Heli_5en.csv', df, 'helix')
    #   df.contact = add_data('BB_Heli_5en.csv',df,'contact')
    fulldata = pd.concat([df_ee, df_rg, df_HB, df_helix])
    fulldata.to_csv('fulldata.csv', index=False)
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


def read_data(directory, filename):
    filepath = os.path.join(directory, filename)
    data = pd.read_csv(filepath)
    return (data)


def select_MTFE_data(data, MTFE):
    MTFEdata = data[data['MTFE'] == MTFE]
    MTFEdata.to_csv('MTFE' + str(MTFE) + 'data.csv', index=False)
    return MTFEdata


def select_datatype_data(data, datatype, MTFE=0):
    typedata = data[data['datatype'] == datatype].reset_index()
    typedata.to_csv('type' + str(datatype) + 'MTFE' + str(MTFE) + 'data.csv', index=False)
    return typedata


def select_protein_data(data, protein):
    protein_data = data[data['Protein'] == protein].reset_index()
    return protein_data

def operation():
    targetmap=[]
    return targetmap
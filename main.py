# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import os

import pandas as pd

import Generatedatabase
import database
import Helical_sliding_window as Hs
def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')
    directory='F:\globus\simulation_contactmap_validation'
    Generatedatabase.generate_database(directory)
    df=pd.read_csv('database_entry.csv')
    # for index, row in df.iterrows():
    #     target=row
        #target=df[df['Protein']=='puma_wildfull']
        # Hs.generate_helicit_csv(target)





def load_data(df):
    df_rg=database.add_data('rg',df)
    df_ee=database.add_data('ee',df)
    fulldata=pd.concat([df_ee,df_rg])
    fulldata.to_csv('fulldata.csv')
    bufferdata=fulldata[fulldata['MTFE']==0.0]
    bufferdata.to_csv('bufferdata.csv')
    eedata=bufferdata[bufferdata['datatype']=='ee'].reset_index()
    eedata.insert(eedata.shape[1],'ratio','Nan')
    rgdata=bufferdata[bufferdata['datatype']=='rg'].reset_index()
    rgdata.to_csv('rgdata.csv')
    for i in range(len(eedata)):
        protein_name=eedata.loc[i,'Protein']
        print(protein_name)
        ee=eedata.loc[i,'Rs']
        rg=rgdata[rgdata['Protein']==protein_name]['Rs']
        print(ee)
        print(rg.iloc[0])
        ratio=ee/rg.iloc[0]
        eedata.loc[i,'ratio']=ratio
    eedata.to_csv('eedata.csv')



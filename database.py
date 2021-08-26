import pandas as pd
import readdata
def add_data(data,entrydf):
    full_df=pd.DataFrame()
    for i in range(len(entrydf)):
        protein_name=entrydf.loc[i,'Protein']
        protein_directory=entrydf.loc[i,'Directory']
        if data=='rg':
            RG=readdata.read_rg(protein_name,protein_directory)
            RG.insert(0,'datatype','rg')
            RG.insert(0,'Protein',protein_name)
            full_df=pd.concat([full_df,RG],ignore_index=True)
        elif data=='ee':
            EE=readdata.read_ee(protein_name,protein_directory)
            EE.insert(0,'datatype','ee')
            EE.insert(0,'Protein',protein_name)
            full_df=pd.concat([full_df,EE],ignore_index=True)
    return full_df

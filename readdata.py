import pandas as pd
import os
def read_rg(protein,directory):
    filename='BBrawRG_5en.csv'
    filepath = os.path.join(directory, filename)
    RG=pd.read_csv(filepath)
    return(RG)
def read_ee(protein,directory):
    filename='BBEERG_5en.csv'
    filepath = os.path.join(directory, filename)
    EE=pd.read_csv(filepath)
    return(EE)
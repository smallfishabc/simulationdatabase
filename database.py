import pandas as pd
import os
import pysnooper

# Load separate datafiles into a entire database
def load_data(df):
    # Add selected data to the dataframe
    # Add Radius of Gyration
    df_rg = add_data('BBrawRG_5en.csv', df, 'rg')
    # Add End to End distance
    df_ee = add_data('BBEERG_5en.csv', df, 'ee')
    # Add average hydrogen bonds number.(The format of HB file is slightly different from other data. This will
    # be fixed in future update)
    df_HB = add_data_HB('BB_HB_5en.csv', df, 'HB')
    # Add Helicity
    df_helix = add_data('BB_Heli_5en.csv', df, 'helix')
    # Add interaction strength
    df_interaction = add_data('BBcontact_new_test1.csv', df, 'interaction')
    # Add sequence features
    df_feature = add_data_sequence_feature(df)
    #df_interaction = add_data('BB_contact_lines.csv',df,'interaction')
    # Merge the data to a single dataframe
    fulldata = pd.concat([df_ee, df_rg, df_HB, df_helix,df_interaction,df_feature])
    # Save the final data to a csv file
    # fulldata.to_csv('fulldata.csv', index=False)
    return fulldata


# Add data function
def add_data(filename, entrydf, data_type):
    # Create an empty dataframe
    full_df = pd.DataFrame()
    # Using loop to add every protein data into the dataframe
    for i in range(len(entrydf)):
        # Retrieve the protein name
        protein_name = entrydf.loc[i, 'Protein']
        # Retrieve the protein directory
        protein_directory = entrydf.loc[i, 'Directory']
        # Read data from the directory
        new_data = read_data(protein_directory, filename)
        # Add datatype to the raw data
        new_data.insert(0, 'datatype', data_type)
        # Add protein name to the raw data
        new_data.insert(0, 'Protein', protein_name)
        # Add these data to the dataframe
        full_df = pd.concat([full_df, new_data], ignore_index=True)
    return full_df

# Because the HB file has different format, we need to rename the HB when attach to the database
# This have a similar structure with previous function and will be removed on next update
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
# This function will use cider and other computational package to calculate the sequence features of the protein
# To be finished
def add_data_sequence_feature(entrydf, data_type='feature'):
    # Create an empty dataframe
    full_df = pd.DataFrame()
    # Using loop to add every protein data into the dataframe
    for i in range(len(entrydf)):
        # Retrieve the protein name
        protein_name = entrydf.loc[i, 'Protein']
        # Retrieve the protein directory
        protein_directory = entrydf.loc[i, 'Directory']
        # Retrieve the protein sequence
        sequence= entrydf.loc[i,'Sequence']
        # Analyse the sequence feature
        length=len(sequence)
        print(sequence,length)
        #Will calculate specific sequence feature
        #cider_calculation()
        #Feng_calculation()
        # Create a new dataframe fo the sequence feature
        new_data=pd.DataFrame({'Protein':protein_name,'datatype':data_type,'Sequence':sequence,'length':length},index=[0])

        #new_data.insert(0, 'datatype', data_type)
        #new_data.insert(0, 'Protein', protein_name)
        # Add these data to the 'full' dataframe
        full_df = pd.concat([full_df, new_data], ignore_index=True)
    return full_df

# Function for retrieving the data from certain directory
def read_data(directory, filename):
    filepath = os.path.join(directory, filename)
    data = pd.read_csv(filepath)
    return (data)



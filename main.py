
import os

import pandas as pd
import Generatedatabase
import database

# If we run this script
if __name__ == '__main__':
    # Set the target directory
    #directory = 'F:\DATA_F\YAP_shahar'
    #directory = 'F:\DATA_F\ADvariants'
    #directory = 'F:\DATA_F\GSlinker_entropic_force'
    #directory = r'F:\DATA_F\PDB_1009'
    #directory = r'F:\DATA_F\LEA_cesar_renamed_2024'
    #directory = 'F:\DATA_F\Entropic_UGDH'
    #directory = 'F:\DATA_F\LEA_ceasar'
    directory = 'F:\DATA_F\Interaction_map_simulation_ACCESS\p53_cancer'
    os.chdir(directory)
    print(directory)
    #print(os.getcwd())
    # Generate a protein entry database containing the subdirectory for each protein
    Generatedatabase.generate_database(directory)
    # Load the entry dataframe
    df = pd.read_csv('database_entry.csv')
    # Load data into the full protein dataframe
    test = database.load_data_easy_no_interation_feature(df)
    #test = database.load_data_easy_entropy(df)
    # Go back to the target directory for saving the csv file
    os.chdir(directory)
    # Save the csv file
    test.to_csv('database_full_value_1119_interaction.csv',index=False)


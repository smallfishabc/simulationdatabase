
import os

import pandas as pd
import correlationplot
import Generatedatabase
import database

# If we run this script
if __name__ == '__main__':
    # Set the target directory
    directory = 'F:\DATA_F\YAP_shahar'
    #directory = 'F:\DATA_F\ADvariants'
    #directory = 'F:\DATA_F\GSlinker_entropic_force'
    #directory = 'F:\DATA_F\LEA_ceasar'
    os.chdir(directory)
    # Special directory setting for linux system
    #linux_directory = '/media/lemoncatboy/WD_BLACK/DATA_F/puma_scramble_new/puma_scrammble_sum'
    #linux_directory = '/media/lemoncatboy/WD_BLACK/DATA_F/PDBsum'
    #directory=linux_directory
    # Print out the directory for testing
    print(directory)
    #print(os.getcwd())
    # Generate a protein entry database containing the subdirectory for each protein
    Generatedatabase.generate_database(directory)
    # Load the entry dataframe
    df = pd.read_csv('database_entry.csv')
    # Load data into the full protein dataframe
    test = database.load_data_easy_no_interation_feature(df)
    #test = database.load_data(df)
    # Go back to the target directory for saving the csv file
    os.chdir(directory)
    # Print for test
    print('a')
    # Save the csv file
    test.to_csv('database_full_value_0428_entropy.csv',index=False)


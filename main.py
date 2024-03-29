
import os

import pandas as pd
import correlationplot
import Generatedatabase
import database
# This simulation database library is used to manage and analyse the large SolSpace intriniscally disorderded Protein (IDP)
# database.

# In this library, each Protein will be simulated under 3-9 different solution condition. Each solution condition will
# have 3-5 repeats. I created a class to automatically detect the simulation configuration and directory of the dataset.

# This library can also generate a dataframe/csv file to summarize the structural property of IDPs including my novel
# interaction map method.

# This is the main function of the library. If we run this script with specific directory, it will save the required data
# to csv.

# Here the database_entry.csv stores the information of simulation settings and the file location.

# The database_full saved the required analysis result.
if __name__ == '__main__':
    # Set the target directory
    directory = 'F:\DATA_F\PDBsum'
    # Special directory setting for linux system
    #linux_directory = '/media/lemoncatboy/WD_BLACK/DATA_F/YAP_shahar'
    #linux_directory = '/media/lemoncatboy/WD_BLACK/DATA_F/puma_scramble_new/puma_scrammble_sum'
    #linux_directory = '/media/lemoncatboy/WD_BLACK/DATA_F/PDBsum'
    #linux_directory = '/media/lemoncatboy/WD_BLACK/DATA_F/PDBsumreal'
    linux_directory = '/media/lemoncatboy/WD_BLACK/DATA_F/LEA_ceasar'
    directory=linux_directory
    # Print out the directory for testing
    print(directory)
    # Generate a Protein entry database containing the subdirectory for each Protein
    Generatedatabase.generate_database(directory)
    # Load the entry dataframe
    df = pd.read_csv('database_entry.csv')
    # Load data into the full Protein dataframe
    test = database.load_data_easy(df)
    # Go back to the target directory for saving the csv file
    os.chdir(directory)
    # Print for test
    print('a')
    # Save the csv file
    #test.to_csv('database_full_value_1103.csv',index=False)


import mdtraj as md
import os
import statistics as st
import numpy as np
import sum_all_easy_Template
import ast
import pandas as pd

# This script is trying to integrate my existing stand alone analysis scripts.
# These stand alone scripts are written in procedural programming.
# We will pass the solution condition using the database_entry file into the sum_easy stand alone analysis script
# We will also read the file location


def read_energy(value_list_origin):
    value_list=ast.literal_eval(value_list_origin)
    print(value_list)
    new_list =[('S_' + x) for x in value_list]
    print(new_list)
    num_list=[]
    num_list =[ast.literal_eval(x) for x in value_list]
    print(num_list)
    return new_list, num_list


def analyze_easy(home_directory, entrydf):
    for i in range(len(entrydf)):
        # Retrieve the protein name
        protein_name = entrydf.loc[i, 'Protein']
        # Retrieve the protein directory
        protein_directory = entrydf.loc[i, 'Directory']
        # Read MTFE from the entry dataframe
        protein_energy, protein_energy_number = read_energy(entrydf.loc[i, 'Psivalue'])
        # Read repeat number from the entry dataframe
        protein_repeat = entrydf.loc[i, 'Repeats']
        # Calculate the ensemble property
        print(protein_energy, protein_energy_number)
        #sum_all_easy_Template.easy_standard(protein_energy, protein_energy_number, protein_repeat, protein_directory)
        # Change back to home directory
        os.chdir(home_directory)

if __name__=="__main__":
    home_directory='F:\DATA_F\GSlinker_entropic_force'
    os.chdir(home_directory)
    df = pd.read_csv('database_entry.csv')
    analyze_easy(home_directory,df)
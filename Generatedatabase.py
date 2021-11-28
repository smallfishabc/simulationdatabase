# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 16:21:42 2021

@author: ShaharGroup-fyu
"""

import os

import pandas as pd

import proteinlistlibrary


# Generate a list of the sub-directory (Old version) (Can be updated next time)
def listsubdirectory(listname, targetdir):
    # Change the directory to the target dir
    os.chdir(targetdir)
    # Generate a full list of folders and files under the directory
    namestemp = os.listdir()
    # Only append folder to the list
    for n in namestemp:
        if (".csv" in n) or (".py" in n) or (".dat" in n) or (".xlsx" in n) or (".json" in n) or (".png" in n) or (
                ".jpg" in n) or ("pycache" in n) or (".svg" in n) or (".txt" in n) or (".pdf" in n) or (".fasta" in n):
            pass
        else:
            listname.append(n)


# Generate the entrance for all Protein trajectory and save their information
def generate_database(string):
    # Confirm the working directory
    print(string)
    # List all sub-directory in the folder
    names = []
    listsubdirectory(names, string)
    # Retrieve the Protein name from the folder name
    disprot_name = []
    for h in names:
        disprot_name.append(h.split("-")[0])  ##payattention no - in folders names
    # Is this line useful?(Double check in next version)
    disprot_name = set(disprot_name)
    # Create a empty list for storing object
    object_list = []
    # Append all Protein object to the list
    for h in disprot_name:
        # Confirm the Protein name
        print(h)
        # Create the Protein object
        q = proteinlistlibrary.Protein(h, string)
        object_list.append(q)
    # Create an empty dataframe
    df = pd.DataFrame(columns=['Protein', 'Directory', 'Sequence', 'Resitype', 'Psivalue', 'Repeats'])
    # Append the data to the empty dataframe
    for i in object_list:
        # To see whether we have single solution condition.
        if len(i.get_sss_solution_type()) == 1:
            # Add required value to the dataframe
            values_to_add = {'Protein': i.name, 'Directory': i.path, 'Sequence': i.get_sequence(),
                             'Resitype': [k.name for k in i.get_sss_solution_type()],
                             'Psivalue': [k.concentration for k in i.type_object],
                             'Repeats': i.type_object[0].repeats}
            row_to_add = pd.DataFrame(values_to_add)
            # Following lines can help me covert this list type to single value in pandas.
            # However, for the convenience, I will still keep the list style here.
            #            row_to_add = row_to_add.explode('Resitype')
            #            row_to_add = row_to_add.explode('Psivalue')
            df = df.append(row_to_add, ignore_index=True)
        else:
            # Currently we only support single solution type, so it will raise a error if we input different conditions.
            raise ValueError('Currently, we only support single solution type')
    # Change back to the target directory and save the csv (Whether necessary? Need to check on next version)
    os.chdir(string)
    df.to_csv('database_entry.csv', index=False)
    return object_list


# This function only generate a object list without generating a dataframe.(Designed to be called by other scripts)
def generate_database_standalone(string):
    # Confirm the working directory
    print(string)
    # List all sub-directory in the folder
    names = []
    listsubdirectory(names, string)
    # Retrieve the Protein name from the folder name
    disprot_name = []
    for h in names:
        disprot_name.append(h.split("-")[0])  ##payattention no - in folders names
    # Is this line useful?(Double check in next version)
    disprot_name = set(disprot_name)
    # Create a empty list for storing object
    object_list = []
    # Append all Protein object to the list
    for h in disprot_name:
        # Confirm the Protein name
        print(h)
        # Create the Protein object
        q = proteinlistlibrary.Protein(h, string)
        object_list.append(q)
    return object_list

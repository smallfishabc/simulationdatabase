# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 16:21:42 2021

@author: ShaharGroup-fyu
"""

import pandas as pd
import numpy as np
import os
import re
import mdtraj as md
import statistics as st
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


# Generate the entrance for all protein trajectory and save their information
def generate_database(string):
    # Confirm the working directory
    print(string)
    # List all sub-directory in the folder
    names = []
    listsubdirectory(names, string)
    # Retrieve the protein name from the folder name
    disprotname = []
    for h in names:
        disprotname.append(h.split("-")[0])  ##payattention no - in folders names
    # Is this line useful?(Double check in next version)
    disprotname = set(disprotname)
    # Create a empty list for storing object
    objectlist = []
    # Append all protein object to the list
    for h in disprotname:
        # Confirm the protein name
        print(h)
        # Create the protein object
        q = proteinlistlibrary.protein(h, string)
        objectlist.append(q)
    # Create an empty dataframe
    df = pd.DataFrame(columns=['Protein', 'Directory', 'Sequence', 'Resitype', 'Psivalue', 'Repeats'])
    # Append the data to the empty dataframe
    for i in objectlist:
        values_to_add = {'Protein': i.name, 'Directory': i.path, 'Sequence': i.getsequence(),
                         'Resitype': [k.name for k in i.getssstype()], 'Psivalue': [k.name for k in i.getssstype()[0].getsssvalue()],
                         'Repeats': i.getssstype()[0].getsssvalue()[0].getrepeats()}
        row_to_add = pd.Series(values_to_add)
        df = df.append(row_to_add, ignore_index=True)
    # Change back to the target directory and save the csv (Whether necessary? Need to check on next version)
    os.chdir(string)
    df.to_csv('database_entry.csv', index=False)
    return objectlist

# This function only generate a object list without generating a dataframe.(Designed to be called by other scripts)
def generate_database_standalone(string):
    # Confirm the working directory
    print(string)
    # List all sub-directory in the folder
    names = []
    listsubdirectory(names, string)
    # Retrieve the protein name from the folder name
    disprotname = []
    for h in names:
        disprotname.append(h.split("-")[0])  ##payattention no - in folders names
    # Is this line useful?(Double check in next version)
    disprotname = set(disprotname)
    # Create a empty list for storing object
    objectlist = []
    # Append all protein object to the list
    for h in disprotname:
        # Confirm the protein name
        print(h)
        # Create the protein object
        q = proteinlistlibrary.protein(h, string)
        objectlist.append(q)
    return objectlist
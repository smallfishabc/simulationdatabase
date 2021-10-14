# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 15:24:30 2021

@author: ShaharGroup-fyu
"""

import pandas as pd
import numpy as np
import os
import re
import mdtraj as md
import statistics as st
import pysnooper


# @pysnooper.snoop()

# This file is a custom tool box for finding subdirectory in each protein folders

# Search file(not folder) in selected directory
def searchfile(root, filename_target):
    # Print out the working directory and target file
    print(root, filename_target)
    # Use os.walk to find the file
    for path, dirs, files in os.walk(root):
        for filename in files:
            if filename == filename_target:
                # Return the full path of the file and jump out of the loop
                fullpath = os.path.join(path, filename)
                break
    return fullpath


# Search folder in selected directory
def search(root, filename):
    # Use os.walk to find the file
    for path, dirs, files in os.walk(root):
        for dirname in dirs:
            if dirname == filename:
                # Return the full path of the file and jump out of the loop
                fullpath = os.path.join(path, dirname)
                break
    return (fullpath)


# Define a class to store the directory information of the protein
class protein:

    # We can add an attribute to actual path
    def __init__(self, name, rootdir=os.getcwd):
        # Set the protein name
        self.name = name
        # Set the protein root directory (Where is the protein folder)
        self.rootdir = rootdir
        # Set the directory of the protein.
        self.path = search(self.rootdir, self.name + '-summary')

    # Get the protein sequence
    def getsequence(self):
        # Double check the protein path
        print(self.path)
        # Search for the seq.fasta file and get its full path
        path = searchfile(self.path, 'seq.fasta')
        # Open the fast file and read the sequence
        seqopen = open(path, 'r')
        seq = seqopen.read()
        # Remove the special character on the tail of the sequence
        while (re.search('\s', seq[-1]) != None):
            seq = seq[:-1]
        return seq

    # Get what type of solution have been simulated for this protein
    def getssstype(self):
        # Create a list for storing types
        typelist = []
        # Find the type by searching the sub-folder
        for path, dirs, files in os.walk(self.path):
            for dirname in dirs:
                print(dirname)
                # Create a type object storing the solution concentration of that type
                dirobject = rtype(dirname, os.path.join(self.path, dirname))
                # Store the type object
                typelist.append(dirobject)
            # add break to avoid further go into the subdirectory (Maxdepth=1)
            break
        return (typelist)


# Define a class for storing the concentration series of the solution
class rtype:
    def __init__(self, name, path):
        # Retrieve the name of the solution
        self.name = name
        # Retrieve the path of the solution folder
        self.path = path

    # Get the concentration information of the solution by searching the sub-folder
    def getsssvalue(self):
        # Create a list for storing concentration
        valuelist = []
        # Find the concentration by searching the sub-folder
        for path, dirs, files in os.walk(self.path):
            for dirname in dirs:
                # Modify the folder name to get the absolute value of the concentration
                value = dirname.split('_', 1)[1]
                # Create a type object storing the number of repeats of that concentration
                dirobject = psssvalue(value, os.path.join(self.path, dirname))
                # Store the concentration object
                valuelist.append(dirobject)
        return (valuelist)

# Define a class for storing the repeats number of certain concentration
class psssvalue:
    def __init__(self, name, path):
        # Retrieve the concentration
        self.name = name
        # Retrieve the path of the folder
        self.path = path
    # Get the repeat number
    def getrepeats(self):
        repeat_number = 0
        # Count the number of xtc file to get the repeat number
        for path, dirs, files in os.walk(self.path):
            for i in files:
                if '.xtc' in i:
                    repeat_number += 1
        return repeat_number

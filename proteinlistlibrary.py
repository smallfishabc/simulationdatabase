# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 15:24:30 2021

@author: ShaharGroup-fyu
"""

import os
import re


# @pysnooper.snoop()

# This file is a custom tool box for finding subdirectory in each Protein folders

# Search file(not folder) in selected directory
def search_file(root, filename_target):
    # Print out the working directory and target file
    print(root, filename_target)
    # Use os.walk to find the file
    for path, dirs, files in os.walk(root):
        for filename in files:
            if filename == filename_target:
                # Return the full path of the file and jump out of the loop
                full_path = os.path.join(path, filename)
                break
    return full_path


# Search folder in selected directory
def search(root, filename):
    # Use os.walk to find the file
    for path, dirs, files in os.walk(root):
        for dirname in dirs:
            if dirname == filename:
                # Return the full path of the file and jump out of the loop
                full_path = os.path.join(path, dirname)
                break
    return full_path


# Define a class to store the directory information of the Protein
class Protein:

    # We can add an attribute to actual path
    def __init__(self, name, rootdir=os.getcwd):
        # Set the Protein name
        self.name = name
        # Set the Protein root directory (Where is the Protein folder)
        self.rootdir = rootdir
        # Set the directory of the Protein.
        self.path = search(self.rootdir, self.name + '-summary')
        # Get the list of solution conditions
        self.type_object = self.get_sss_solution_type()

    # Get the Protein sequence
    def get_sequence(self):
        # Double check the Protein path
        print(self.path)
        # Search for the seq.fasta file and get its full path
        path = search_file(self.path, 'seq.fasta')
        # Open the fast file and read the sequence
        seqopen = open(path, 'r')
        seq = seqopen.read()
        # Remove the special character on the tail of the sequence
        while re.search('\s', seq[-1]) is not None:
            seq = seq[:-1]
        return seq

    # Get what type of solution have been simulated for this Protein
    def get_sss_solution_type(self):
        # Create a list for storing types objects
        object_list = []
        # Find the type by searching the sub-folder
        for path, dirs, files in os.walk(self.path):
            for dirname in dirs:
                # Create a type object storing the solution concentration of that type
                type_object = SolutionType(dirname, os.path.join(self.path, dirname))
                # Store the type object
                object_list.append(type_object)
            # add break to avoid further go into the subdirectory (Max_depth=1)
            break
        print([k.name for k in object_list])
        return object_list


# Define a class for storing the concentration series of the solution
class SolutionType:
    def __init__(self, name, path):
        # Retrieve the name of the solution
        self.name = name
        # Retrieve the path of the solution folder
        self.path = path
        # Retrieve the concentration list of the solution folder
        self.concentration, self.location = self.get_concentration_value()
        # Retrieve the repeats number of this condition (Here we assume the repeats number is the same across all
        # concentrations)
        self.repeats = self.get_repeats()

    # Get the concentration information of the solution by searching the sub-folder
    def get_concentration_value(self):
        # Create a list for storing concentration
        value_list = []
        # Create a list for storing the path of concentration folders
        path_list = []
        # Find the concentration by searching the sub-folder
        for path, dirs, files in os.walk(self.path):
            for dirname in dirs:
                # Modify the folder name to get the absolute value of the concentration
                value = dirname.split('_', 1)[1]
                # Return the full path of the concentration directory
                dir_full_path = os.path.join(self.path, dirname)
                # Store the concentration value
                value_list.append(value)
                # Store the path of the concentration folder
                path_list.append(dir_full_path)
        return value_list, path_list

    # Get the repeat number
    def get_repeats(self):
        repeat_number = 0
        # Count the number of xtc file to get the repeat number
        for path, dirs, files in os.walk(self.location[0]):
            for i in files:
                if '.xtc' in i:
                    repeat_number += 1
        return repeat_number

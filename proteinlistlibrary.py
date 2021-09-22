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

#@pysnooper.snoop()
def searchfile(root, filename):
    print(root, filename)
    for path, dirs, files in os.walk(root):
        for dirname in files:
            if dirname == filename:
                fullpath = os.path.join(path, dirname)
                break
    return (fullpath)


def search(root, filename):
    for path, dirs, files in os.walk(root):
        for dirname in dirs:
            if dirname == filename:
                fullpath = os.path.join(path, dirname)
                break
    return (fullpath)


class protein:
    # We can add an attribute to actual path
    def __init__(self, name, rootdir=os.getcwd):
        self.name = name
        self.rootdir = rootdir
        self.path = search(self.rootdir, self.name + '-summary')

    def getssstype(self):
        typelist = []
        for path, dirs, files in os.walk(self.path):
            for dirname in dirs:
                print(dirname)
                dirobject = rtype(dirname, os.path.join(self.path, dirname))
                typelist.append(dirobject)
            #add break to avoid further go into the subdirectory
            break
        return (typelist)

    def getsequence(self):
        print(self.path)
        path = searchfile(self.path, 'seq.fasta')
        seqopen = open(path, 'r')
        seq = seqopen.read()
        while (re.search('\s', seq[-1]) != None):
            seq = seq[:-1]
        return (seq)


class rtype:
    def __init__(self, name, path):
        self.name = name
        self.path = path

    def getsssvalue(self):
        valuelist = []
        for path, dirs, files in os.walk(self.path):
            for dirname in dirs:
                value = dirname.split('_', 1)[1]
                dirobject = psssvalue(value, os.path.join(self.path, dirname))
                valuelist.append(dirobject)
        return (valuelist)


class psssvalue:
    def __init__(self, name, path):
        self.name = name
        self.path = path

    def getrepeats(self):
        repeatnumber = 0
        for path, dirs, files in os.walk(self.path):
            for i in files:
                if '.xtc' in i:
                    repeatnumber += 1
        return (repeatnumber)

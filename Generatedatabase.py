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
def listsubdirectory(listname,targetdir):   
    os.chdir(targetdir)
    namestemp = os.listdir()
    for n in namestemp:
        if (".csv"in n) or (".py" in n) or (".dat" in n) or (".xlsx" in n) or (".json" in n) or(".png" in n) or (".jpg" in n) or ("pycache" in n) or (".svg" in n):
            pass
        else:
            listname.append(n)
    os.chdir(string)
names=[]
string=str(os.getcwd())
listsubdirectory(names,string)
disprotname=[]
for h in names:
    disprotname.append(h.split("-")[0]) ##payattention no - in folders names
disprotname=set(disprotname)
objectlist=[]
for h in disprotname:
    q= proteinlistlibrary.protein(h, string)
    objectlist.append(q)
df=pd.DataFrame(columns=['Protein','Directory','Sequence','Resitype','Psivalue','Repeats'])
for i in objectlist:
    values_to_add = {'Protein': i.name,'Directory': i.path,'Sequence': i.getsequence(),'Resitype': i.getssstype(),'Psivalue': i.getssstype()[0].getsssvalue(),'Repeats':i.getssstype()[0].getsssvalue()[0].getrepeats()}
    row_to_add = pd.Series(values_to_add)
    df = df.append(row_to_add,ignore_index=True)
df.to_csv('database_entry.csv', index=False)
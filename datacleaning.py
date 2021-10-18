import os

import pandas as pd
import ast
def covert_string_tolist(stringlist):
    string_cleaned = ast.literal_eval(stringlist)
    return(string_cleaned)

def absolute_value(list):
    new_list=[]
    for i in list:
        absvalue=abs(i)
        new_list.append(absvalue)
    return new_list
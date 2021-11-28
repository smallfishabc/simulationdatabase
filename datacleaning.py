import ast
# Convert the list in csv file back to a list in python
def covert_string_tolist(stringlist):
    string_cleaned = ast.literal_eval(stringlist)
    return(string_cleaned)
# Get the absolute value of the list
def absolute_value(list):
    new_list=[]
    for i in list:
        abs_value=abs(i)
        new_list.append(abs_value)
    return new_list
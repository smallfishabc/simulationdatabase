# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import os

import pandas as pd
import correlationplot
import Generatedatabase
import database


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')
    directory = 'F:\DATA_F\PDBsum'
    Generatedatabase.generate_database(directory)
    df = pd.read_csv('database_entry.csv')
    test = database.load_data(df)
    os.chdir(directory)
    print('a')
    test.to_csv('database_full_value_test2.csv',index=False)
#    operation = database.operation(df)
#    correlationplot.custom_plot(df)


import os

import pandas as pd

import datacleaning


# Here, we will define several database operation that are frequently used in my previous scripts.

# Select concentration we want
def select_mtfe_data(data, mtfe):
    mtfe_data = data[data['MTFE'] == mtfe]
    return mtfe_data


# Select the data type we want (Rg,Re,...) Default concentration will be 0 (water)
def select_datatype_data(data, datatype):
    type_data = data[data['datatype'] == datatype].reset_index()
    return type_data


# Select the Protein we want
def select_protein_data(data, protein):
    protein_data = data[data['Protein'] == protein].reset_index()
    return protein_data


# Make several selection at the same time (Will modify this according to my scripts)
def multiple_selection_function(data, protein, datatype, mtfe):
    # First we select by the protein name
    df1 = select_protein_data(data, protein)
    # Then we select by the datatype(Rg,Ree...)
    df2 = select_datatype_data(df1, datatype)
    # Then we select by the MTFE value
    df3 = select_mtfe_data(df2, mtfe)
    return df3


# Calculate the Re of ideal polymer
def GS_length_calculation(sequence):
    length = len(sequence)
    # These pre-factors are the fitting result of GS_linker simulation
    return 0.55359657 * (length) ** 0.47499067


# Calculate the Re Ratio between target Protein and ideal polymer (Normalized Ree called Chi)
def chi_value_calculation(data, protein, mtfe, sequence):
    ee = multiple_selection_function(data, protein, 'ee', mtfe)['Rs'].tolist()[0]
    l_gs = GS_length_calculation(sequence)
    chi = ee / l_gs - 1
    return chi


# Calculate the solution sensitivity of proteins in attractive solutions(Chi method)
def attractive_sensitivity_calculation(data, protein, sequence, chi_0):
    chi_p3 = chi_value_calculation(data, protein, 3, sequence)
    sensitivity = chi_p3 - chi_0
    return sensitivity


# Calculate the solution sensitivity of proteins in repulsive solutions(Chi method)
def repulsive_sensitivity_calculation(data, protein, sequence, chi_0):
    chi_m3 = chi_value_calculation(data, protein, -3, sequence)
    sensitivity = chi_m3 - chi_0
    return sensitivity


# Calculate the solution sensitivity of proteins in attractive solutions(Ree method)
def attractive_sensitivity_calculation_ree(data, protein, datatype='ee'):
    ee_p3 = multiple_selection_function(data, protein, datatype, 3)
    ee_0 = multiple_selection_function(data, protein, datatype, 0)
    sensitivity = ee_p3 - ee_0
    return sensitivity


# Calculate the solution sensitivity of proteins in repulsive solutions(Ree method)
def repulsive_sensitivity_calculation_ree(data, protein, datatype='ee'):
    ee_m3 = multiple_selection_function(data, protein, datatype, -3)
    ee_0 = multiple_selection_function(data, protein, datatype, 0)
    sensitivity = ee_m3 - ee_0
    return sensitivity


# Designed for calculating several value in on function
def database_plot_pre_multiple(datatype, mtfe=0, index=8):
    # os.chdir('/media/lemoncatboy/WD_BLACK/DATA_F/puma_scramble_new/puma_scrammble_sum')
    # Read the location of simulation file
    df_entry = pd.read_csv('database_entry.csv')
    # Read the structural information of the protein
    df = pd.read_csv('database_full_value_1103.csv')
    # Create an empty Dataframe for storing the information
    full_df = pd.DataFrame()
    # Looping over every protein
    for i in range(len(df_entry)):
        # Get protein's name
        protein_name = df_entry.loc[i, 'Protein']
        # Get portein's sequence
        sequence = df_entry.loc[i, 'Sequence']
        # Add sequence length for normalization (If we do not need normalization, we put 1 here)
        seq_length = 1
        # Add the df_protein for selecting sequence feature
        df_protein = df[df['Protein'] == protein_name]
        # Covert list from csv to python list and get attractive and repulsive interaction
        # These interaction are defined in my interaction map project
        valueraw1_raw = multiple_selection_function(df, protein_name, datatype, mtfe)['att1'].tolist()[0]
        valueraw1_list = datacleaning.covert_string_tolist(valueraw1_raw)
        valueraw2_raw = multiple_selection_function(df, protein_name, datatype, mtfe)['att2'].tolist()[0]
        valueraw2_list = datacleaning.covert_string_tolist(valueraw2_raw)
        valueraw3_raw = multiple_selection_function(df, protein_name, datatype, mtfe)['rep1'].tolist()[0]
        valueraw3_list = datacleaning.covert_string_tolist(valueraw3_raw)
        # Here the repulsive interaction are in minus value so we need to get their absolute value
        valueraw3_list = datacleaning.absolute_value(valueraw3_list)
        valueraw4_raw = multiple_selection_function(df, protein_name, datatype, mtfe)['rep2'].tolist()[0]
        valueraw4_list = datacleaning.covert_string_tolist(valueraw4_raw)
        valueraw4_list = datacleaning.absolute_value(valueraw4_list)
        # Here, we are trying different interaction calculation algorithm.
        # We stored all possible interaction types in the typelist
        # Record the datatype of each value using typelist
        type_list = ['none_none', 'none_ratio', 'none_probability', 'dis_none', 'dis_ratio', 'dis_probability',
                    'rtdis_none', 'rtdis_ratio', 'rtdis_probability']

        # Fast calculation of interaction strength using single index. Can be removed in later version
        # type_inlist = type_list[index]
        # valueraw1 = valueraw1_list[index]
        # valueraw2 = valueraw2_list[index]
        # valueraw3 = valueraw3_list[index]
        # valueraw4 = valueraw4_list[index]
        # value = (valueraw1 +  valueraw2 + valueraw3 + valueraw4)
        # value = (valueraw1 + valueraw2 - valueraw3 - valueraw4)
        # value_att = (valueraw1 + valueraw2)
        # value_rep = (valueraw3 + valueraw4)
        # print(value)

        # Pass these interaction list to other function to generate a panda data frame
        interaction_data = interaction_pre_list(index, type_list, seq_length, valueraw1_list, valueraw2_list,
                                                valueraw3_list, valueraw4_list)
        # Calculate the chi value of this protein
        chi_0 = chi_value_calculation(df, protein_name, 0, sequence)
        # Calculate the solution sensitivity of the protein
        att_sensitivity = attractive_sensitivity_calculation(df, protein_name, sequence, chi_0)
        rep_sensitivity = repulsive_sensitivity_calculation(df, protein_name, sequence, chi_0)
        # Calculate the helicity ot the protein
        heli = multiple_selection_function(df, protein_name, 'helix', 0)['Rs'].tolist()[0]
        # Calculate the H-bonds fraction of the protein
        hb = multiple_selection_function(df, protein_name, 'HB', 0)['Rs'].tolist()[0]
        # Calculate the beta-sheet fraction of the protein
        beta = multiple_selection_function(df, protein_name, 'beta', 0)['Rs'].tolist()[0]
        # Calculate Kappa value
        delta = select_datatype_data(df_protein, 'feature', mtfe=0)['delta']
        # Create a diction containing all these structural information
        dict_pd = {'Protein': protein_name, 'chi_0': chi_0, 'att_sensitivity': att_sensitivity,
                   'rep_sensitivity': rep_sensitivity, 'Helicity': heli, 'H-bonds': hb, 'delta': delta, 'beta': beta}
        # Covert the dictionary to the dataframe
        new_data = pd.DataFrame(dict_pd, index=[0])
        # Combine the structural information with the interaction data
        new_data = pd.concat([new_data, interaction_data], axis=1)
        # Combine the existing data with this protein's data
        full_df = pd.concat([full_df, new_data], ignore_index=True)
    # Save the final result to csv file
    full_df.to_csv('interaction_strength_fitting_1014.csv', index=False)
    return full_df

# Here we covert the interaction type list to a full pandas data frame
def plot_pre_list(datatype):
    list_data = []
    # How many different interaction type should be taken into consideration
    for i in range(9):
        j = database_plot_pre_multiple(datatype, mtfe=0, index=i)
        list_data.append(j)
    return list_data

# Here we covert the interaction type list to a full pandas data frame
def interaction_pre_list(index, typelist, seq_length, valueraw1_list, valueraw2_list, valueraw3_list, valueraw4_list):
    full_data = pd.DataFrame({'a': 1}, index=[0])
    for i in range(index + 1):
        # Retrieve value from the list
        valueraw1 = valueraw1_list[i]
        valueraw2 = valueraw2_list[i]
        valueraw3 = valueraw3_list[i]
        valueraw4 = valueraw4_list[i]
        # Calculate the normalized interaction strength
        rawatt = (valueraw1 + valueraw2) / seq_length
        rawrep = (valueraw3 + valueraw4) / seq_length
        # Add column name to the dataframe
        column_name = typelist[i]
        column_name_att = column_name + 'att'
        column_name_rep = column_name + 'rep'
        full_data[column_name_att] = rawatt
        full_data[column_name_rep] = rawrep
    full_data.drop(columns=['a'], inplace=True)
    return full_data


#####################
# Here is the backup Code for debugging
# A plot function designed for finding correlation between different datatype.
def database_plot_pre(datatype, mtfe=0):
    df_entry = pd.read_csv('database_entry.csv')
    df = pd.read_csv('database_full_value_test1.csv')
    full_df = pd.DataFrame()
    for i in range(len(df_entry)):
        protein_name = df_entry.loc[i, 'Protein']
        sequence = df_entry.loc[i, 'Sequence']
        valueraw1 = multiple_selection_function(df, protein_name, datatype, mtfe)['att1'].tolist()[0]
        valueraw2 = multiple_selection_function(df, protein_name, datatype, mtfe)['att2'].tolist()[0]
        valueraw3 = multiple_selection_function(df, protein_name, datatype, mtfe)['rep1'].tolist()[0]
        valueraw4 = multiple_selection_function(df, protein_name, datatype, mtfe)['rep2'].tolist()[0]
        # These lines provides different possible interaction strength.
        # value=(valueraw1+2*valueraw2-valueraw3-2*valueraw4)/len(sequence)
        # value = (valueraw1 + 2 * valueraw2) / len(sequence)
        # value = (valueraw1 + valueraw2) / len(sequence)
        value = (valueraw1 + valueraw2 + valueraw3 + valueraw4) / len(sequence)
        chi_0 = chi_value_calculation(df, protein_name, 0, sequence)
        att_sensitivity = attractive_sensitivity_calculation(df, protein_name, sequence, chi_0)
        rep_sensitivity = repulsive_sensitivity_calculation(df, protein_name, sequence, chi_0)
        dict_pd = {'Protein': protein_name, datatype: value, 'chi_0': chi_0, 'att_sensitivity': att_sensitivity,
                   'rep_sensitivity': rep_sensitivity}
        new_data = pd.DataFrame(dict_pd, index=[0])
        full_df = pd.concat([full_df, new_data], ignore_index=True)
    return full_df


#########
# Here the main function is designed for debugging
if __name__ == '__main__':
    os.chdir('/media/lemoncatboy/WD_BLACK/DATA_F/PDBsum')
    database_plot_pre_multiple('interaction', mtfe=0)

import mdtraj as md
import os
import statistics as st
import numpy as np
import sum_all_easy_Template
import entropy_easy_0502_template
import entropy_easy_0502_enhanced_sampling_template
import entropy_easy_0517_sphere_template
import entropy_easy_0517_sphere_reverse_template
import platdist_0504_template
import ast
import pandas as pd

# This script is trying to integrate my existing stand alone analysis scripts.
# These stand alone scripts are written in procedural programming.
# We will pass the solution condition using the database_entry file into the sum_easy stand alone analysis script
# We will also read the file location


def read_energy(value_list_origin):
    value_list=ast.literal_eval(value_list_origin)
    print(value_list)
    new_list =[('S_' + x) for x in value_list]
    print(new_list)
    num_list=[]
    num_list =[ast.literal_eval(x) for x in value_list]
    print(num_list)
    return new_list, num_list


def analyze_easy(home_directory, entrydf):
    for i in range(len(entrydf)):
        # Retrieve the protein name
        protein_name = entrydf.loc[i, 'Protein']
        # Retrieve the protein directory
        protein_directory = entrydf.loc[i, 'Directory']
        # Read MTFE from the entry dataframe
        protein_energy, protein_energy_number = read_energy(entrydf.loc[i, 'Psivalue'])
        # Read repeat number from the entry dataframe
        protein_repeat = entrydf.loc[i, 'Repeats']
        # Calculate the ensemble property
        print(protein_energy, protein_energy_number)
        os.chdir(protein_directory)
        #if os.path.exists('BBforbid_d0.2_0504.csv'):
            #os.rename('BBforbid_d0.2_0504.csv','BBforbid_d0_0504_backup.csv')
        #if os.path.exists('BBrawRG_easy_0430.csv'):
        #    print(i)
        #    continue
        #else:
            #sum_all_easy_Template.easy_standard(protein_energy, protein_energy_number, protein_repeat, protein_directory)
        #distance_angle_pair_list=[(0,np.pi/6),(0.2,np.pi/6),(0.5,np.pi/6),(1,np.pi/6)]
        #radius_list=[1,3,10,20]
        radius_list=[0, 0.2, 0.5, 1, 2]
        #radius_list = [4 ,5, 6]
        #radius_list = [7 ,8, 9]
        #radius_list = [12,14,16,18]
        #radius_list = [11,13,17,19]
        #radius_list = [15,20,25,30]
        #radius_list = [35, 40, 45, 50]
        #radius_list = [75, 100, 150, 200]
        # for i in distance_angle_pair_list:
        for i in radius_list:
            if os.path.exists('BB_entropy_easy_'+'radius'+str(i)+'_0517_reverse.csv'):
                print(i)
                continue
            else:
                entropy_easy_0517_sphere_reverse_template.sphere_entropy_reverse(protein_energy,protein_energy_number,protein_repeat, protein_directory,i)
            # if os.path.exists('BB_entropy_easy_'+'radius'+str(i)+'_0517.csv'):
            #     print(i)
            #     continue
            # else:
            #     entropy_easy_0517_sphere_template.sphere_entropy(protein_energy,protein_energy_number,protein_repeat, protein_directory,i)
            # entropy_easy_0502_template.cone_entropy(protein_energy, protein_energy_number, protein_repeat,
            #                                         protein_directory, i[0], i[1])
            #if os.path.exists('BB_entropy_enhanced_easy_'+'d'+str(round(i[0],2))+'_theta'+str(round(i[1],2))+'_0502.csv'):
                #print(i)
                #continue
            #else:
                #entropy_easy_0502_enhanced_sampling_template.enhanced_sampling_surface(protein_energy,protein_energy_number,protein_repeat, protein_directory,i[0], i[1])
            # if os.path.exists('BBforbid_d'+str(round(i[0],2))+'_0504.csv'):
            #     print(i)
            #     continue
            # else:
            #     platdist_0504_template.plat_entropy(protein_energy,protein_energy_number,protein_repeat, protein_directory,i[0])
        # Change back to home directory
        os.chdir(home_directory)

if __name__=="__main__":
    #home_directory='F:\DATA_F\PDBsumreal'
    home_directory='F:\DATA_F\GSlinker_entropic_force'
    os.chdir(home_directory)
    df = pd.read_csv('database_entry.csv')
    analyze_easy(home_directory,df)
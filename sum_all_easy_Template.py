#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 16:24:47 2021

@author: lemoncatboy
"""
import mdtraj as md
import os, sys
import statistics as st
import numpy as np
import pandas as pd
import re

# Calculate the Hydrogen-bond number using Mdtraj library
# Here the r is the trajectory that only containing the protein atoms
# j is the number of frames in the trajectory.
def calc_HB(r,j,length):
    #Calculate the Hydogen_bond
    hbo=md.wernet_nilsson(r)
    # This is a couter for the loop over frames
    op=0
    # A list to store the final result.
    HB=[]
    # looping over each frame
    while op<j:
        # hbo[op] is a sub numpy array stored donor and acceptor 
        # information of the H-bonds. The length of the array is the number of
        # H-bonds
        if (len(hbo[op])):
            HB.append(len(hbo[op])/length)
        op+=1
    # Average H-bonds per residue per frame
    return(sum(HB)/j)

# Calculate the helical propensity fo the protein.
# Here the t is the full trajectory file(We can try to use r instead)
def calc_Heli(t,length):
    # Calculate the secondary structure using dssp
    dssp = md.compute_dssp(t)
    # Create an empty array with the row number equals to the number of residues
    # We may not need the number 1 here
    dssp_count = np.zeros((1, t.n_residues))
    # Loop over frames
    for i in range(t.n_frames):
        # Loop over residues
        for j in range(t.n_residues):
            # In dssp, letter 'H' means alpha-helix
            if dssp[i,j] in 'H':
                dssp_count[0,j] += 1
    # Get the average of the dssp_count
    dssp_prob = np.divide(dssp_count,t.n_frames)
    summary=dssp_prob.sum(axis=1)
    return(summary[0]/length)

#calc_beta is added 03/28/2022
# Same algorithm as alpha
def calc_Beta(t,length):
    dssp = md.compute_dssp(t)
    dssp_count = np.zeros((1, t.n_residues))
    for i in range(t.n_frames):
        for j in range(t.n_residues):
            if dssp[i,j] in 'E':
                dssp_count[0,j] += 1
    dssp_prob = np.divide(dssp_count,t.n_frames)
    summary = dssp_prob.sum(axis=1)
    return(summary[0]/length)

# Average radius of gyration of the ensemble
def calc_Rg(r):
    d = md.compute_rg(r)
    return(st.mean(d))
# Average end to end distance of the ensemble
# The defination of the end to end distnace is the distance between 
# the first alpha carbon and the last alpha carbon of the protein sequence.

def calc_Re(r,t):
    topology=t.topology
    rpology=topology.select_atom_indices(selection='alpha')
    d = md.compute_distances(r,[[rpology[0],rpology[-1]]])
    listtemp=[]
    for temp in d:
        listtemp.append(float(temp[0]))
    return(st.mean(listtemp))
 
# Get the average of five trajectories. 
# Here, we calculate the weighted average and standard deviation by
# frame numbers     
def average_function(mean_value,stlist,frametraj_temp,meantraj_temp,number_of_repeats):
    framesum=0
    meanadj=0
    for nn in range(0,number_of_repeats,1):
        framesum+=frametraj_temp[nn]
        meanadj+=frametraj_temp[nn]*meantraj_temp[nn]
    meanref=meanadj/framesum
    mean_value.append(meanref)
    stdadj=0
    for nn in range(0,number_of_repeats,1):
        stdadj+=frametraj_temp[nn]*(meantraj_temp[nn]-meanref)*(meantraj_temp[nn]-meanref)
    stlist.append(np.sqrt((stdadj/framesum)*number_of_repeats/(number_of_repeats-1)))

#Read sequence from the fast file
def read_seq():
    seqopen=open('seq.fasta','r')
    seq=seqopen.read()
    while (re.search('\s',seq[-1]) != None):
        seq=seq[:-1]
    x=len(seq)
    return seq,x

# A protocol function for calculating all quantities above.
def easy_standard(k,q,repeat,pwd):
    os.chdir(pwd)
    print (pwd)
    seq,length=read_seq()
    h = 'BB'
    mean_HB=[]
    mean_Heli=[]
    mean_Beta=[]
    mean_Rg=[]
    mean_Re=[]
    st_HB=[]
    st_Heli=[]
    st_Beta=[]
    st_Rg=[]
    st_Re=[]
    frametraj=[]
    for p in k:
        meantraj_temp=[[] for i in range(repeat)]
        frametraj_temp=[]
        string = str(pwd)+'/'+h+'/'+p
        os.chdir(string)
        test = os.getcwd()
        print (test)    
        for nn in range (0,repeat,1):
            try:
                t = md.load({'__traj_'+str(nn)+'.xtc'},top='__START_0.pdb')
            except:
                t = md.load({'__traj_'+str(nn)+'.xtc'},top='__END_0.pdb')
            u=t.top.select('protein')
            r=t.atom_slice(u)
            j = r.n_frames
            frametraj_temp.append(j)
            
            HB_temp=calc_HB(r,j,length)
            Heli_temp=calc_Heli(t,length)
            Beta_temp=calc_Beta(t,length)
            Rg_temp=calc_Rg(r)
            Re_temp=calc_Re(r,t)
            
            meantraj_temp[0].append(HB_temp)
            meantraj_temp[1].append(Heli_temp)
            meantraj_temp[2].append(Beta_temp)
            meantraj_temp[3].append(Rg_temp)
            meantraj_temp[4].append(Re_temp)
            
        
        average_function(mean_HB,st_HB,frametraj_temp,meantraj_temp[0],repeat)
        average_function(mean_Heli,st_Heli,frametraj_temp,meantraj_temp[1],repeat)
        average_function(mean_Beta,st_Beta,frametraj_temp,meantraj_temp[2],repeat)
        average_function(mean_Rg,st_Rg,frametraj_temp,meantraj_temp[3],repeat)
        average_function(mean_Re,st_Re,frametraj_temp,meantraj_temp[4],repeat)
        frametraj.append(sum(frametraj_temp))
        
    os.chdir(pwd)    
    dataframe0 = pd.DataFrame({'MTFE':q,'Hbond':mean_HB,'st':st_HB,'Frame':frametraj})
    dataframe0.sort_values(by=['MTFE'],inplace=True)
    pcsv0=h+'_HB_easy_0430.csv'
    dataframe0.to_csv(pcsv0,index=False,sep=',')
    dataframe1 = pd.DataFrame({'MTFE':q,'Rs':mean_Heli,'Sd':st_Heli,'Frame':frametraj})
    dataframe1.sort_values(by=['MTFE'],inplace=True)
    pcsv1=h+'_Heli_easy_0430.csv'
    dataframe1.to_csv(pcsv1,index=False,sep=',')
    dataframe2 = pd.DataFrame({'MTFE':q,'Rs':mean_Beta,'Sd':st_Beta,'Frame':frametraj})
    dataframe2.sort_values(by=['MTFE'],inplace=True)
    pcsv2=h+'_Beta_easy_0430.csv'
    dataframe2.to_csv(pcsv2,index=False,sep=',')
    dataframe3 = pd.DataFrame({'MTFE':q,'Rs':mean_Rg,'Sd':st_Rg,'Frame':frametraj})
    dataframe3.sort_values(by=['MTFE'],inplace=True)
    pcsv3=h+'rawRG_easy_0430.csv'
    dataframe3.to_csv(pcsv3,index=False,sep=',')
    dataframe4 = pd.DataFrame({'MTFE':q,'Rs':mean_Re,'Sd':st_Re,'Frame':frametraj})
    dataframe4.sort_values(by=['MTFE'],inplace=True)
    pcsv4=h+'EERG_easy_0430.csv'
    dataframe4.to_csv(pcsv4,index=False,sep=',')

if __name__=="__main__":
    k = ['S_-3','S_-2','S_-1','S_0','S_1','S_2','S_3']
    q =[-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0]
    repeat=5
    pwd=os.path.dirname(os.path.realpath(__file__))
    easy_standard(k,q,repeat,pwd)

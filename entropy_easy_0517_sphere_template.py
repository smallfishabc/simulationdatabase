#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 16:24:47 2021

@author: lemoncatboy
"""

import mdtraj as md
import os,sys
import numpy as np
import pandas as pd
import re
import entropy_library as el
import shape_protocol_library as spl

#Read sequence from the fast file
def read_seq():
    seqopen=open('seq.fasta','r')
    seq=seqopen.read()
    while (re.search('\s',seq[-1]) != None):
        seq=seq[:-1]
    x=len(seq)
    return seq,x

# A protocol function for calculating entropy with a sphere constraint.
def sphere_entropy(k,q,repeat,pwd,radius):
    os.chdir(pwd)
    print (pwd)
    seq,length=read_seq()
    h = 'BB'
    
    omega_list=[]
    frametraj=[]
    
    for p in k:
        string = str(pwd)+'/'+h+'/'+p
        os.chdir(string)
        test = os.getcwd()
        print (test)
        traj_list=['__traj_'+str(nn)+'.xtc' for nn in range (0,repeat,1)]
        try:
            t = md.load(traj_list,top='__START_0.pdb')
        except:
            t = md.load(traj_list,top='__END_0.pdb')
        u=t.top.select('protein')
        r=t.atom_slice(u)
        j = r.n_frames
        topology=t.topology
        r_alpha=topology.select_atom_indices(selection='alpha')
        prohibited,ratio=spl.compute_forbidden_curvature(radius,t,r_alpha,j)
        prohibited_df = pd.DataFrame(prohibited)
        prohibited_df.to_csv('prohibited'+'radius'+str(radius)+'_frame.csv', index=False, sep=',')
        # Record forbidden ratio and frame number
        omega_list.append(ratio/j)
        frametraj.append(j)
    os.chdir(pwd)
    dataframe = pd.DataFrame({'MTFE':q,'Omega2/Omega1':omega_list,'Frame':frametraj})
    pcsv=h+'_entropy_easy_'+'radius'+str(radius)+'_0517.csv'
    dataframe.sort_values(by=['MTFE'], inplace=True)
    dataframe.to_csv(pcsv,index=False,sep=',')


if __name__=="__main__":
    k = ['S_-3','S_-2','S_-1','S_0','S_1','S_2','S_3']
    q =[-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0]
    repeat=5
    #pwd=os.path.dirname(os.path.realpath(__file__))
    pwd ='F:\DATA_F\GSlinker_entropic_force\GS8-summary'
    radius=1
    sphere_entropy_entropy(k,q,repeat,pwd,radius)

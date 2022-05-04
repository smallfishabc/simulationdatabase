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

# A protocol function for calculating entropy with a surface constratint 
# and enhanced sampling algorithm.
# Compared to other scripts, here the angle_theta = (np.pi-angle_theta)
def enhanced_sampling_surface(k,q,repeat,pwd,distance_d,angle_theta):
    os.chdir(pwd)
    print (pwd)
    seq,length=read_seq()
    h = 'BB'
    ratio_list=[]
    frametraj=[]
    
    for p in k:
        string = str(pwd)+'/'+h+'/'+p
        os.chdir(string)
        test = os.getcwd()
        print(test)
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
        # The angle defination is different so we need to use np.pi/2-angle_theta
        result=spl.compute_forbiden_rotation_enhanced(distance_d,np.pi/2-angle_theta,t,r_alpha,j)
        ratio=np.average(result)
        # Record forbidden ratio and frame number
        ratio_list.append(ratio)
        frametraj.append(j)
    os.chdir(pwd)    
    dataframe = pd.DataFrame({'MTFE':q,'Omega2/Omega1':ratio_list,'Frame':frametraj})
    pcsv=h+'_entropy_enhanced_easy_'+'d'+str(round(distance_d,2))+'_theta'+str(round(angle_theta,2))+'_0502.csv'
    dataframe.sort_values(by=['MTFE'], inplace=True)
    dataframe.to_csv(pcsv,index=False,sep=',')

if __name__=="__main__":
    k = ['S_-3','S_-2','S_-1','S_0','S_1','S_2','S_3']
    q =[-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0]
    repeat=5
    pwd=os.path.dirname(os.path.realpath(__file__))
    distance_d=0
    angle_theta=0.001
    enhanced_sampling_surface(k,q,repeat,pwd,distance_d,angle_theta)

# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 13:22:07 2019

@author: ShaharGroup-fyu
"""

import matplotlib.pyplot as plt
import mdtraj as md
import os, sys
import statistics as st
import numpy as np
import pandas as pd
import re
def read_seq():
    seqopen=open('seq.fasta','r')
    seq=seqopen.read()
    while (re.search('\s',seq[-1]) != None):
        seq=seq[:-1]
    x=len(seq)
    return seq,x

def computeforbiden(distance_d,t,r_alpha,jframes):
    k=0
    j=0
    prohibited_frames=[]
    while j<jframes:
        for i in r_alpha:
            valuedot=1
            value=1
            standard=t.xyz[j,r_alpha[0],:]
            standardvector=t.xyz[j,r_alpha[1],:]-t.xyz[j,r_alpha[0],:]
            a=t.xyz[j,i,:]
            if i>r_alpha[1]:
                #valuea=(standardvector[0]-standard[0])*(standardvector[0]-standard[0])+(standardvector[1]-standard[1])*(standardvector[1]-standard[1])+(standardvector[2]-standard[2])*(standardvector[2]-standard[2])
                valuea=np.linalg.norm(standardvector)
                #valueb=a-standard
                #valueb=(a[0]-standard[0])*(a[0]-standard[0])+(a[1]-standard[1])*(a[1]-standard[1])+(a[2]-standard[2])*(a[2]-standard[2])
                #valuedot=(standardvector[0]-standard[0])*(a[0]-standard[0])+(standardvector[1]-standard[1])*(a[1]-standard[1])+(standardvector[2]-standard[2])*(a[2]-standard[2])
                valuedot=np.dot(standardvector,(a-standard))
                value=valuedot/np.sqrt(valuea)
            if value<=distance_d:
                prohibited_frames.append(j)
                k+=1
                break
        j+=1   
    return(1-k/jframes,prohibited_frames)
arrange=0
col=0
def plat_entropy(k,q,repeat,pwd,distance_d):
    os.chdir(pwd)
    print (pwd)
    seq,length=read_seq()
    h = 'BB'
    omega_list=[]
    frametraj=[]
    for p in k:
        string = str(pwd)+'/'+h+'/'+p
        os.chdir(string)
        print(string)
        traj_list=['__traj_'+str(nn)+'.xtc' for nn in range (0,repeat,1)]
        try:
            t = md.load(traj_list,top='__START_0.pdb')
        except:
            t = md.load(traj_list,top='__END_0.pdb')
        j = t.n_frames
        topology=t.topology
        r_alpha=topology.select_atom_indices(selection='alpha')
        frametraj.append(j)
        #-1 means lower than the platform
        distance_d_adjust=-1*distance_d
        temp,prohibited=computeforbiden(distance_d_adjust,t,r_alpha,j)
        prohibited_df = pd.DataFrame(prohibited)
        prohibited_df.to_csv('prohibited'+str(distance_d_adjust)+'_frame.csv', index=False, sep=',')
        omega_list.append(temp)
    os.chdir(pwd)
    dataframe = pd.DataFrame({'MTFE':q,'Omega2/Omega1':omega_list,'Frame':frametraj})
    pcsv=h+'forbid_d'+str(distance_d)+'_0504.csv'
    dataframe.sort_values(by=['MTFE'], inplace=True)
    dataframe.to_csv(pcsv,index=False,sep=',')

if __name__=="__main__":
    k = ['S_-3','S_-2','S_-1','S_0','S_1','S_2','S_3']
    q =[-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0]
    repeat=5
    #pwd=os.path.dirname(os.path.realpath(__file__))
    pwd='F:\DATA_F\GSlinker_entropic_force\GS16-summary'
    distance_d=1
    plat_entropy(k,q,repeat,pwd,distance_d)        

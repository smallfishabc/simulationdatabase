import os
import statistics

import mdtraj as md
import numpy as np
import matplotlib.pylab as plt
import calculation
import pandas as pd

# Sliding window script for future data analysis. May be removed in next update

def helicity_calculation(repeats, directory):
    os.chdir(directory)
    test=os.getcwd()
    print(test)
    frametraj = []
    meantraj = []
    for nn in range(repeats):
        t = md.load({'__traj_' + str(nn) + '.xtc'}, top='__START_0.pdb')
        u = t.top.select('Protein')
        r = t.atom_slice(u)
        j = r.n_frames
        frametraj.append(j)
        d = md.compute_dssp(r)
        translate=translate_helicity(d)
        meantraj.append(translate)
    mean,sd,frame=calculation.mean_feng(meantraj,frametraj,5)
    return (mean, sd, frame)


def translate_helicity(d):
    output = np.zeros(d.shape)
    #illiterate array(better)
    for index1,i in enumerate(d):
        for index2,j in enumerate(i):
            if j == 'H':
                output[index1][index2] = 1
            else:
                output[index1][index2] = 0
    average=output.mean(axis=0)
    #
    return (average)


def sliding_window(mean):
    #Generator (better)
    size=5
    returnlist=[]
    returnvaluelist=[]
    chunklist=[]
    chunkvaluelist=[]
    for index,i in enumerate(mean):
        if index>size/2-1 and index<len(mean)-size/2-1:
            chunk=[index-2,index-1,index,index+1,index+2]
            chunkvalue=[mean[index-2],mean[index-1],i,mean[index+1],mean[index+2]]
            chunklist.append(chunk)
            chunkvaluelist.append(chunkvalue)
            returnlist.append(index)
            returnvaluelist.append(statistics.mean(chunkvalue))
        i+=1
    return returnlist,returnvaluelist



def cidertype_plot(x,y,title):
    plt.scatter(x, y, s=50, marker='o', color='Black', zorder=5)
    plt.xlabel('Helical propensity(blob5)')
    plt.ylabel('Residue index')
    plt.title(title)
    plt.ylim(0,1)
    plt.show()

def generate_helicit_csv(target):
    mean,sd,frame=helicity_calculation(5,target['Directory']+'\BB\S_0')
    sequence=target['Sequence']
    x,y=sliding_window(mean)
    cidertype_plot(x,y,'PUMA_scramble_2')
    a=[None]*len(sequence)
    b=list(range(1,len(sequence)+1))
    for index,i in enumerate(x):
        a[i]=y[index]
    for index , i in enumerate(a):
        if i==None:
            a[index]=0
    os.chdir(target['Directory'])
    savedf=pd.DataFrame({'residue index': b,'amino acid': sequence,'helicity(blob 5)': a})
    savedf.to_csv(target['Protein']+'helicity0910.csv',index=None,header='None')

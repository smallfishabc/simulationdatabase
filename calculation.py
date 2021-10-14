import mdtraj as md
import os
import statistics as st
import numpy as np

# This script is trying to integrate my existing stand alone analysis scripts.
# These stand alone scripts are written in procedural programming.

# Calculate the Radius of Gyration of the trajectory
def calculate_rg(repeats,directory):
    os.chdir(directory)
    frametraj=[]
    meantraj=[]
    for nn in range(repeats):
        t = md.load({'__traj_' + str(nn) + '.xtc'}, top='__START_0.pdb')
        u = t.top.select('protein')
        r = t.atom_slice(u)
        j = r.n_frames
        frametraj.append(j)
        d = md.compute_rg(r)
        meantraj.append(st.mean(d))
    mean,sd,frame=mean_feng(meantraj,frametraj)
    return (mean,sd,frame)

# Different repeats may have different number of frames. We create a algorithm to calculated the average based on its
# frame number
def mean_feng(value,count,repeats):
    framesum=0
    meanadj=0
    for nn in range(repeats):
        framesum += count[nn]
        meanadj += count[nn] * value[nn]
    meanref = meanadj / framesum
    stdadj = 0
    for nn in range(repeats):
        stdadj += count[nn] * (value[nn] - meanref) * (value[nn] - meanref)
    sd=np.sqrt((stdadj / framesum) * repeats / (repeats-1))
    return (meanref,sd,framesum)

# Calculation template for these stand alone functions.
def calculation_structure(protein,function):
    string=protein.path
    os.chdir(string)
    l=protein.getssstype()
    for h in l:
        mean = []
        sd = []
        frame = []
        k=h.getsssvalue()
        for p in k:
            repeats=p.getrepeats()
            runpath = str(string) + '/' + h.name + '/' + p.name
            if function=='rg':
                meanvalue,sdvalue,framevalue=calculate_rg(repeats,runpath)
            mean.append(meanvalue)
            sd.append(sdvalue)
            frame.append(framevalue)
    return()
# Here, we will add a function to automatically run the existing scripts under the certain folder.
def run_exist_scripts(datatype,entry_directory):
    os.chdir(entry_directory)
    if datatype is 'abc' :
        pass
    elif 1:
        pass
    elif 1:
        pass
    elif 1:
        pass
    return()

def chi_calculation():
    return

def std_calculation():
    return

def temp_calculation():
    return
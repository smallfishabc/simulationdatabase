import mdtraj as md
import os
import statistics as st
import numpy as np

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


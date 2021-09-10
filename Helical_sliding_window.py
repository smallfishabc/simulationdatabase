import os

import mdtraj as md
import numpy as np

import calculation


def helicity_calculation(repeats, directory):
    os.chdir(directory)
    test=os.getcwd()
    print(test)
    frametraj = []
    meantraj = []
    for nn in range(repeats):
        t = md.load({'__traj_' + str(nn) + '.xtc'}, top='__START_0.pdb')
        u = t.top.select('protein')
        r = t.atom_slice(u)
        j = r.n_frames
        frametraj.append(j)
        d = md.compute_dssp(r)
        meantraj.append(translate_helicity(d))
        mean,sd,frame=calculation.mean_feng(meantraj,frametraj)
    return (mean, sd, frame)


def translate_helicity(d):
    output = np.zeros(d.shape)
    #illiterate array(better)
    for index1,j in enumerate(i):
        for index2,k in enumerate(j):
            if k is 'H':
                output[index1][index2] = 1
            else:
                output[index1][index2] = 0
    average=output.mean(axis=1)
    #
    return (average)


def sliding_window():
    return


def cidertype_plot():
    return

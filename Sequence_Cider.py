import pandas as pd
import numpy as np
from localcider.sequenceParameters import SequenceParameters

# Create a Cider sequence object for sequence parameter calculation
def create_Seq_Object(sequence):
    SeqOb = SequenceParameters(sequence)
    return SeqOb
# We need to pass a Seq_object to calculate following sequence feature
# Calculate the kappa value which describe the charge residue dispersion
def calculate_kappa(SeqOb):
    kappa=SeqOb.get_kappa()
    return kappa
# Calculate the fraction of charge residue in the sequence
def calculate_FCR(SeqOb):
    FCR=SeqOb.get_FCR()
    return FCR
# Calculate the net charge per residue in the sequence (+ for positive - for negative)
def calculate_NCPR(SeqOb):
    NCPR=SeqOb.get_NCPR()
    return NCPR
# Calculate the hydrophobicity of the sequence
def calculate_mean_hydropathy(SeqOb):
    mean_hydropathy=SeqOb.get_mean_hydropathy()
    return mean_hydropathy
# Calculate the fraction of expanding residue(E/D/R/K/P) in the sequence
def calculate_fraction_expanding(SeqOb):
    fraction_expanding=SeqOb.get_fraction_expanding()
    return fraction_expanding
# Calculate the charge residue dispersion raw value
def calculate_delta(SeqOb):
    delta=SeqOb.get_delta()
    return delta

# This function is to calculate all parameters and return to main function with a dictionary
def Cider_calculation(sequence):
    SeqOb=create_Seq_Object(sequence)
    kappa=calculate_kappa(SeqOb)
    FCR=calculate_FCR(SeqOb)
    NCPR=calculate_NCPR(SeqOb)
    mean_hydropathy=calculate_mean_hydropathy(SeqOb)
    fraction_expanding=calculate_fraction_expanding(SeqOb)
    delta=calculate_delta(SeqOb)
    dict_parameter={'kappa':kappa,'FCR':FCR,'NCPR':NCPR,'Hydrophobicity':mean_hydropathy
                    ,'Expanding':fraction_expanding,'delta':delta}
    return dict_parameter
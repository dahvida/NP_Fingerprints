from map4 import MAP4Calculator
from mhfp.encoder import MHFPEncoder
import numpy as np

######################################################################

MHFP = MHFPEncoder(1024)
MAP4 = MAP4Calculator(dimensions=1024) 

######################################################################

def jaccard_like(a,b):
    return float(np.count_nonzero(a == b)) / float(len(a))

def calc_MAP4(mols):
    fps = MAP4.calculate_many(mols)
    return np.array(fps)

def calc_MHFP(mols):
    fps = [MHFP.encode_mol(x) for x in mols]
    return np.array(fps)



"""MinHashing Fingerprint functions

Used in fp_script.py
Based on:
https://github.com/reymond-group/mhfpk
https://github.com/reymond-group/map4
"""

from FP_calc.map4 import MAP4Calculator
from FP_calc.mhfp import MHFPEncoder
import numpy as np
from typing import *

############################################################################

#initialize calculators with default FP size
MHFP = MHFPEncoder(1024)
MAP4 = MAP4Calculator(dimensions=1024) 

############################################################################

def efficient_array(fp_function):
    """
    Decorator for all MinHashing functions, so that they take as input
    rdkit mols and return numpy arrays
    """
    def wrapper(*args):
        print(f"[FP]: Executing {fp_function.__name__}")
        raw_fps = fp_function(*args)

        #have to store in double precision due to large integers in FPs
        #from MinHashing
        array = np.array(raw_fps, dtype=np.int32)

        return array
    
    return wrapper

#############################################################################

@efficient_array
def calc_MAP4(
        mols: List
        ) -> List:
    return MAP4.calculate_many(mols)

#--------------------------------------------------------------------------#

@efficient_array
def calc_MHFP(
        mols: List
        ) -> List:
    return [MHFP.encode_mol(x) for x in mols]



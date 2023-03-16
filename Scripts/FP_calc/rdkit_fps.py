"""RDKIT Fingerprint functions

Used in fp_script.py
Based on: https://www.rdkit.org/docs/
"""

import rdkit
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys
import numpy as np
from rdkit.Chem.AtomPairs import Pairs, Torsions
from rdkit.Avalon.pyAvalonTools import GetAvalonCountFP
from typing import *

############################################################################

def efficient_array(fp_function):
    """
    Decorator for all RDKIT Fingerprint functions, so that they take as input
    rdkit mols and return numpy arrays, with efficient conversion from
    sparse bit vectors into arrays.
    """
    def wrapper(*args):
        print(f"[FP]: Executing {fp_function.__name__}")

        #get fps, number of molecules and FP size
        fps, n_mols, n_bits = fp_function(*args)

        #create empty array of correct size
        array = np.empty((n_mols, n_bits), dtype=np.int16)
        
        #fill array with FPs
        for i in range(n_mols):
            DataStructs.ConvertToNumpyArray(fps[i], array[i])
        
        return array
    
    return wrapper

############################################################################

@efficient_array
def calc_MACCS(
        mols: List[rdkit.Chem.rdchem.Mol]
        ) -> Tuple[List, int, int]:
    fps = [MACCSkeys.GenMACCSKeys(x) for x in mols]
    n_mols = len(fps)
    n_bits = len(fps[0])
    return fps, n_mols, n_bits

#--------------------------------------------------------------------------#

@efficient_array
def calc_ECFP(
        mols: List[rdkit.Chem.rdchem.Mol],
        radius: int = 2,
        nbits: int = 1024
        ) -> Tuple[List, int, int]:
    fps = [AllChem.GetMorganFingerprintAsBitVect(
        x, radius, nbits) for x in mols]
    n_mols = len(fps)
    return fps, n_mols, nbits

#--------------------------------------------------------------------------#

@efficient_array
def calc_FCFP(
        mols: List[rdkit.Chem.rdchem.Mol],
        radius: int = 2,
        nbits: int = 1024
        ) -> Tuple[List, int, int]:
    fps = [AllChem.GetMorganFingerprintAsBitVect(
        x, radius, nbits, useFeatures=True) for x in mols]
    n_mols = len(fps)
    return fps, n_mols, nbits

#--------------------------------------------------------------------------#

@efficient_array
def calc_RDKIT(
        mols: List[rdkit.Chem.rdchem.Mol]
        ) -> Tuple[List, int, int]:
    fps = [Chem.RDKFingerprint(x) for x in mols]
    n_mols = len(fps)
    n_bits = len(fps[0])
    return fps, n_mols, n_bits

#--------------------------------------------------------------------------#

@efficient_array
def calc_AP(
        mols: List[rdkit.Chem.rdchem.Mol],
        nbits: int = 2048
        ) -> Tuple[List, int, int]:
    fps = [Pairs.GetHashedAtomPairFingerprint(x, nbits) for x in mols]
    n_mols = len(fps)
    return fps, n_mols, nbits

#--------------------------------------------------------------------------#

@efficient_array
def calc_TT(
        mols: List[rdkit.Chem.rdchem.Mol],
        nbits: int = 2048
        ) -> Tuple[List, int, int]:
    fps = [Torsions.GetHashedTopologicalTorsionFingerprint(
        x, nbits) for x in mols]
    n_mols = len(fps)
    return fps, n_mols, nbits

#--------------------------------------------------------------------------#

@efficient_array
def calc_AVALON(
        mols: List[rdkit.Chem.rdchem.Mol],
        nbits: int = 512
        ) -> Tuple[List, int, int]:
    fps = [GetAvalonCountFP(x, nbits) for x in mols]
    n_mols = len(fps)
    return fps, n_mols, nbits




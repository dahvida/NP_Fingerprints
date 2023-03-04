import rdkit
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys
import numpy as np
from rdkit.Chem.AtomPairs import Pairs, Torsions
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
from rdkit.Avalon.pyAvalonTools import GetAvalonCountFP

def efficient_array(fp_function):
    def wrapper(*args):
        fps, n_mols, n_bits = fp_function(*args)
        array = np.empty((n_mols, n_bits))
        for i in range(n_mols):
            DataStructs.ConvertToNumpyArray(fps[i], array[i])
        return array
    return wrapper

@efficient_array
def calc_MACCS(mols):
    fps = [MACCSkeys.GenMACCSKeys(x) for x in mols]
    n_mols = len(fps)
    n_bits = len(fps[0])
    return fps, n_mols, n_bits

@efficient_array
def calc_ECFP(mols, radius = 2, nbits = 1024):
    fps = [AllChem.GetMorganFingerprintAsBitVect(x, radius, nbits) for x in mols]
    n_mols = len(fps)
    return fps, n_mols, nbits

@efficient_array
def calc_FCFP(mols, radius = 2, nbits = 1024):
    fps = [AllChem.GetMorganFingerprintAsBitVect(x, radius, nbits, useFeatures = True) for x in mols]
    n_mols = len(fps)
    return fps, n_mols, nbits

@efficient_array
def calc_RDKIT(mols):
    fps = [Chem.RDKFingerprint(x) for x in mols]
    n_mols = len(fps)
    n_bits = len(fps[0])
    return fps, n_mols, n_bits

@efficient_array
def calc_AP(mols, nbits = 2048):
    fps = [Pairs.GetHashedAtomPairFingerprint(x, nbits) for x in mols]
    n_mols = len(fps)
    return fps, n_mols, nbits

@efficient_array
def calc_TT(mols, nbits = 2048):
    fps = [Torsions.GetHashedTopologicalTorsionFingerprint(x, nbits) for x in mols]
    n_mols = len(fps)
    return fps, n_mols, nbits

@efficient_array
def calc_Avalon(mols, nbits = 512):
    fps = [GetAvalonCountFP(x, nbits) for x in mols]
    n_mols = len(fps)
    return fps, n_mols, nbits































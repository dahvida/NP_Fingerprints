import os
import numpy as np
from jpype import java, isJVMStarted, startJVM, getDefaultJVMPath, JPackage

"""
Based on: https://github.com/hcji/pycdk
"""

######################################################################

startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=%s" % "./cdk-2.2.jar")
cdk = JPackage('org').openscience.cdk

def parse_smiles(smiles):
    function = cdk.smiles.SmilesParser(cdk.DefaultChemObjectBuilder.getInstance())
    smiles = [function.parseSmiles(x) for x in smiles] 
    return smiles

def unpack_fps(fps):
    nbits = fps[0].size()
    n_mols = len(fps)
    array = np.zeros((n_mols, nbits), dtype=np.float16)
    
    for i in range(n_mols):
        fp = fps[i].asBitSet()
        bits = []
        idx = fp.nextSetBit(0)
        while idx >= 0:
            bits.append(idx)
            idx = fp.nextSetBit(idx + 1)
        bits = np.array(bits)
        array[i, bits] = 1
    return array

def efficient_array(fp_function):
    def wrapper(*args):
        parsed_smiles = parse_smiles(*args)
        raw_fps = fp_function(parsed_smiles)
        array = unpack_fps(raw_fps)
        return array
    return wrapper

######################################################################

@efficient_array
def calc_PUBCHEM(smiles):
    function_fp = cdk.fingerprint.PubchemFingerprinter(cdk.silent.SilentChemObjectBuilder.getInstance())
    fps = [function_fp.getBitFingerprint(x) for x in smiles]
    return fps
    
@efficient_array
def calc_CDK(smiles, size=1024, depth=7):
    function_fp = cdk.fingerprint.Fingerprinter(size, depth)
    fps = [function_fp.getBitFingerprint(x) for x in smiles]
    return fps

@efficient_array
def calc_ESTATE(smiles):
    function_fp = cdk.fingerprint.EStateFingerprinter()
    fps = [function_fp.getBitFingerprint(x) for x in smiles]
    return fps

@efficient_array
def calc_KR(smiles):
    function_fp = cdk.fingerprint.KlekotaRothFingerprinter()
    fps = [function_fp.getBitFingerprint(x) for x in smiles]
    return fps

@efficient_array
def calc_LINGO(smiles):
    function_fp = cdk.fingerprint.LingoFingerprinter()
    fps = [function_fp.getBitFingerprint(x) for x in smiles]
    return fps


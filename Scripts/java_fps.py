import os
import numpy as np
from jpype import java, isJVMStarted, startJVM, getDefaultJVMPath, JPackage

"""
Based on: https://github.com/hcji/pycdk
"""

startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=%s" % "./cdk-2.2.jar")
cdk = JPackage('org').openscience.cdk

def parse_smiles(smiles):
    function = cdk.smiles.SmilesParser(cdk.DefaultChemObjectBuilder.getInstance())
    smiles = [function.parseSmiles(x) for x in smiles] 
    return smiles

def unpack_fps(fps):
    nbits = fps[0].size()
    n_mols = len(fps)
    array = np.zeros((n_mols, nbits))
    
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

def calc_PUBCHEM(smiles):
    function_fp = cdk.fingerprint.PubchemFingerprinter(cdk.silent.SilentChemObjectBuilder.getInstance())
    fps = [function_fp.getBitFingerprint(x) for x in smiles]
    return fps
    
def calc_CDK(smiles, size, depth):
    function_fp = cdk.fingerprint.Fingerprinter(size, depth)
    fps = [function_fp.getBitFingerprint(x) for x in smiles]
    return fps

def calc_ESTATE(smiles):
    function_fp = cdk.fingerprint.EStateFingerprinter()
    fps = [function_fp.getBitFingerprint(x) for x in smiles]
    return fps

def calc_KR(smiles):
    function_fp = cdk.fingerprint.KlekotaRothFingerprinter()
    fps = [function_fp.getBitFingerprint(x) for x in smiles]
    return fps









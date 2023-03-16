"""CDK Fingerprint functions

Used in fp_script.py
Based on: https://github.com/hcji/pycdk
"""

import os
import numpy as np
from jpype import java, isJVMStarted, startJVM, getDefaultJVMPath, JPackage
from typing import *

############################################################################

#start JAVA instance in python
startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=%s" % "./FP_calc/cdk-2.2.jar")
cdk = JPackage('org').openscience.cdk

#--------------------------------------------------------------------------#

def parse_smiles(
        smiles: List[str],
        addH: bool = False
        ) -> List[cdk.AtomContainer]:
    """Parses SMILES into CDK molecules

    Args:
        smiles: SMILES (M,) to parse    
        addH:   whether to add mols or not
    
    Returns:
        List (M,) of CDK objects
    """
    #create CDK parser instance
    function = cdk.smiles.SmilesParser(cdk.DefaultChemObjectBuilder.getInstance())

    #parse all SMILES in list as IAtomContainer objects
    smiles = [function.parseSmiles(x) for x in smiles]

    #optionally add explicit hydrogen to mols
    if addH is True:

        #create converter
        addH = cdk.tools.manipulator.AtomContainerManipulator.convertImplicitToExplicitHydrogens

        #loop over all mols
        for smile in smiles:
            addH(smile)
    
    return smiles

#--------------------------------------------------------------------------#

def unpack_fps(
        fps: List
        ) -> np.ndarray:
    """Converts CDK fingerprint in numpy array
    
    Args:
        fps:    CDK fingerprints (M,K) to convert

    Returns:
        Fingerprint as numpy array (M,K)
    """
    
    #get FP size K
    nbits = fps[0].size()

    #get number of mols M
    n_mols = len(fps)

    #construct array where FPs will be stored
    array = np.zeros((n_mols, nbits), dtype=np.int16)
    
    #loop over every mol
    for i in range(n_mols):
        
        #convert raw FP into set of all present bits
        fp = fps[i].asBitSet()

        #initialize first bit position
        bits = []

        #store first bit
        idx = fp.nextSetBit(0)

        #append all bits for i-th FP
        while idx >= 0:
            bits.append(idx)
            idx = fp.nextSetBit(idx + 1)
        bits = np.array(bits)
        
        #Check for failures (happened with KR)
        try:
            #store FP into array
            array[i, bits] = 1
        except:
            pass

    return array

#--------------------------------------------------------------------------#

def efficient_array(fp_function):
    """
    Decorator for all CDK FP functions, so that they take as input
    a list of SMILES and output a numpy array.
    """
    def wrapper(*args):
        print(f"[FP]: Executing {fp_function.__name__}")
        
        #convert SMILES in CDK objects
        parsed_smiles = parse_smiles(*args)
        
        #get FPs from CDK objects
        raw_fps = fp_function(parsed_smiles)
        
        #store FPs in numpy arrays
        array = unpack_fps(raw_fps)
        
        return array
    
    return wrapper

############################################################################

@efficient_array
def calc_PUBCHEM(
        smiles: List,
        addH: bool = True
        ) -> List:
    #According to documentation, needs explicit hydrogens
    #https://cdk.github.io/cdk/2.3/docs/api/org/openscience/cdk/fingerprint/PubchemFingerprinter.html
    function_fp = cdk.fingerprint.PubchemFingerprinter(cdk.silent.SilentChemObjectBuilder.getInstance())
    fps = [function_fp.getBitFingerprint(x) for x in smiles]
    return fps

#--------------------------------------------------------------------------#    

@efficient_array
def calc_DAYLIGHT(
        smiles: List,
        addH: bool = True,
        size: int = 1024,
        depth: int = 7
        ) -> List:
    #According to documentation, needs explicit hydrogens
    #http://cdk.github.io/cdk/2.2/docs/api/org/openscience/cdk/fingerprint/Fingerprinter.html
    function_fp = cdk.fingerprint.Fingerprinter(size, depth)
    fps = [function_fp.getBitFingerprint(x) for x in smiles]
    return fps

#--------------------------------------------------------------------------#

@efficient_array
def calc_ESTATE(
        smiles: List
        ) -> List:
    function_fp = cdk.fingerprint.EStateFingerprinter()
    fps = [function_fp.getBitFingerprint(x) for x in smiles]
    return fps

#--------------------------------------------------------------------------#

@efficient_array
def calc_KR(
        smiles: List
        ) -> List:
    function_fp = cdk.fingerprint.KlekotaRothFingerprinter()
    fps = [function_fp.getBitFingerprint(x) for x in smiles]
    return fps

#--------------------------------------------------------------------------#

@efficient_array
def calc_LINGO(
        smiles: List
        ) -> List:
    function_fp = cdk.fingerprint.LingoFingerprinter()
    fps = [function_fp.getBitFingerprint(x) for x in smiles]
    return fps


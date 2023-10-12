"""JCompoundMapper Fingerprint functions

Used in fp_script.py
Based on: https://github.com/hcji/pycdk
"""

import os
import numpy as np
from jpype import java, isJVMStarted, startJVM, getDefaultJVMPath, JPackage
from typing import *
import numpy as np

############################################################################

startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=%s" % "./FP_calc/jCMapperCLI.jar")
jmap = JPackage('de').zbit.jcmapper.fingerprinters
cdk = JPackage('org').openscience.cdk

#--------------------------------------------------------------------------#

def parse_smiles(
        smiles: List[str]
        ) -> List[cdk.Molecule]:
    """Parses SMILES into CDK molecules

    Args:
        smiles: SMILES (M,) to parse        
    Returns:
        List (M,) of CDK objects
    """
    #create CDK parser instance
    function = cdk.smiles.SmilesParser(cdk.DefaultChemObjectBuilder.getInstance())

    #parse all SMILES in list as IAtomContainer objects
    smiles = [function.parseSmiles(x) for x in smiles]
        
    return smiles

#--------------------------------------------------------------------------#

def unpack_fps(
        fps: List,
        ) -> np.ndarray:
    """Converts CDK fingerprint in numpy array
    
    Args:
        fps:    CDK fingerprints (M,K) to convert

    Returns:
        Fingerprint as numpy array (M,K)
    """
    
    #get FP size K
    nbits = 4096

    #get number of mols M
    n_mols = len(fps)

    #construct array where FPs will be stored
    array = np.zeros((n_mols, nbits), dtype=np.int16)
    
    #loop over every mol
    for i in range(n_mols):

        #initialize first bit position
        bits = []

        #store first bit
        idx = fps[i].nextSetBit(0)

        #append all bits for i-th FP
        while idx >= 0:
            bits.append(idx)
            idx = fps[i].nextSetBit(idx + 1)
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
def calc_DFS(
        smiles: List
        ) -> List:
    function_fp = jmap.topological.DepthFirstSearch()
    fps = [function_fp.getHashFingerprint(x,7,4096) for x in smiles]
    
    return fps

@efficient_array
def calc_ASP(
        smiles: List
        ) -> List:
    function_fp = jmap.topological.Encoding2DAllShortestPath()
    fps = [function_fp.getFingerprint(x) for x in smiles]
    fps = [jmap.features.FeatureMap(x).getHashedFingerPrint(4096) for x in fps]

    return fps

@efficient_array
def calc_LSTAR(
        smiles: List
        ) -> List:
    function_fp = jmap.topological.Encoding2DLocalAtomEnvironment()
    fps = [function_fp.getFingerprint(x) for x in smiles]
    fps = [jmap.features.FeatureMap(x).getHashedFingerPrint(4096) for x in fps]

    return fps

@efficient_array
def calc_RAD2D(
        smiles: List
        ) -> List:
    function_fp = jmap.topological.Encoding2DMolprint()
    fps = [function_fp.getFingerprint(x) for x in smiles]
    fps = [jmap.features.FeatureMap(x).getHashedFingerPrint(4096) for x in fps]

    return fps

@efficient_array
def calc_PH2(
        smiles: List
        ) -> List:
    function_fp = jmap.topological.Encoding2DPharmacophore2Point()
    fps = [function_fp.getFingerprint(x) for x in smiles]
    fps = [jmap.features.FeatureMap(x).getHashedFingerPrint(4096) for x in fps]

    return fps

@efficient_array
def calc_PH3(
        smiles: List
        ) -> List:
    function_fp = jmap.topological.Encoding2DPharmacophore3Point()
    fps = [function_fp.getFingerprint(x) for x in smiles]
    fps = [jmap.features.FeatureMap(x).getHashedFingerPrint(4096) for x in fps]

    return fps


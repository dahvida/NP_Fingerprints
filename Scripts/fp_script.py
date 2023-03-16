"""Fingerprint calculation script

Script to compute fingerprints for the preprocessed COCONUT database.
All calculated fingerprint arrays will be pickled into ../FPs as
'fp_codename.pkl', i.e. 'ecfp.pkl'.

Before running this script, make sure you have executed cleanup_script.py

Due to jpype's behavior, it is impossible to import simultaneously
CDK and JCompoundMapper fingerprint functions. As such, the user has
to specify which fingerprint package he wants to use for the calculations
via the --FP_type argument.
"""

from utils import pickle_save
from rdkit import Chem
import pandas as pd
from FP_calc.rdkit_fps import *
from FP_calc.minhash_fps import *
import argparse

###########################################################################

parser = argparse.ArgumentParser(description=__doc__,
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--FP_type",
                    default = "rdkit",
                    help = "fingerprinting package to use for calculations, options are [rdkit, cdk, minhash, jmap]")

args = parser.parse_args()

###########################################################################

def main(fp_class):
    
    fp_class = fp_class.lower()
    print("[FP]: Starting fingerprint generation")
    print(f"[FP]: Will compute fingerprints from {fp_class} package and store in ../FPs")

    #load dataset, get smiles and mol objects
    db = pd.read_csv("../Data/clean_coconut.csv")
    smiles = list(db["sugar_free_smiles"])
    mols = [Chem.MolFromSmiles(x) for x in smiles]
    print("[FP]: Dataset loaded")
    
    #run and save rdkit FPs
    if fp_class == "rdkit":
        fp = calc_ECFP(mols)
        pickle_save(fp, "../FPs/ecfp.pkl", verbose=False)

        fp = calc_AVALON(mols)
        pickle_save(fp, "../FPs/avalon.pkl", verbose=False)

        fp = calc_TT(mols)
        pickle_save(fp, "../FPs/tt.pkl", verbose=False)
 
        fp = calc_MACCS(mols)
        pickle_save(fp, "../FPs/maccs.pkl", verbose=False)

        fp = calc_FCFP(mols)
        pickle_save(fp, "../FPs/fcfp.pkl", verbose=False)

        fp = calc_AP(mols)
        pickle_save(fp, "../FPs/ap.pkl", verbose=False)

        fp = calc_RDKIT(mols)
        pickle_save(fp, "../FPs/rdkit.pkl", verbose=False)
    
    #run and save minhash FPs
    elif fp_class == "minhash":
        fp = calc_MAP4(mols)
        pickle_save(fp, "../FPs/map4.pkl", verbose=False)

        fp = calc_MHFP(mols)
        pickle_save(fp, "../FPs/mhfp.pkl", verbose=False)

    #run and save cdk FPs
    elif fp_class == "cdk":
        from FP_calc.cdk_fps import calc_PUBCHEM, calc_DAYLIGHT, calc_KR, calc_LINGO, calc_ESTATE

        fp = calc_PUBCHEM(smiles)
        pickle_save(fp, "../FPs/pubchem.pkl", verbose=False)

        fp = calc_DAYLIGHT(smiles)
        pickle_save(fp, "../FPs/daylight.pkl", verbose=False)

        fp = calc_KR(smiles)
        pickle_save(fp, "../FPs/kr.pkl", verbose=False)

        fp = calc_LINGO(smiles)
        pickle_save(fp, "../FPs/lingo.pkl", verbose=False)

        fp = calc_ESTATE(smiles)
        pickle_save(fp, "../FPs/estate.pkl", verbose=False)

    #run and save jmap FPs
    elif fp_class == "jmap":
        from FP_calc.jmap_fps import calc_DFS, calc_ASP, calc_LSTAR, calc_RAD2D, calc_PH2, calc_PH3

        fp = calc_DFS(smiles)
        pickle_save(fp, "../FPs/dfs.pkl", verbose=False)

        fp = calc_ASP(smiles)
        pickle_save(fp, "../FPs/asp.pkl", verbose=False)

        fp = calc_LSTAR(smiles)
        pickle_save(fp, "../FPs/lstar.pkl", verbose=False)
        
        fp = calc_RAD2D(smiles)
        pickle_save(fp, "../FPs/rad2d.pkl", verbose=False)
        
        fp = calc_PH2(smiles)
        pickle_save(fp, "../FPs/ph2.pkl", verbose=False)

        fp = calc_PH3(smiles)
        pickle_save(fp, "../FPs/ph3.pkl", verbose=False)
    
    #catch wrong input from user
    else:
        print("[FP]: FP type not recognized, try one of [rdkit, minhash, cdk, jmap]")

    print("[FP]: Calculation finished")



if __name__ == "__main__":
    main(args.FP_type)








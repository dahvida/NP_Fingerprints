"""Utility functions for cleanup_script.py

Preprocessing functions to clean the raw COCONUT database.
SMILES standardization and salt removal is handled by the ChEMBL pipeline package.
label_text is roughly adapted from https://github.com/reymond-group/Coconut-TMAP-SVM
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
from chembl_structure_pipeline import standardizer
from typing import *

###########################################################################

def label_text(
        string: str
        ) -> Tuple[str, int]:
    """Converts taxonomy string into labels

    Args:
        string:     raw taxonomy string from COCONUT_DB

    Returns:
        A tuple containing (1) the organism of origin and (2) the organism
        numerical identifier
    """
    #put everything in lowercase for consistency when matching
    string = string.lower()
    
    #load labelling presets from Coconut-TMAP-SVM
    bacteria_check = ["bacteria", "bacillus", "bacta", "bacterium"]
    plants_check = ["plants", "plant"]
    fungi_check = ["fungi", "fungus", "aspergillus"]
    animal_check = ["animal", "animals"]
    marine_check = ["marine"]
    
    idx_list = []
    
    #run checks and assign label
    #(some NPs have multiple sources but we are assigning according
    #to the first label)
    if any([x in string for x in bacteria_check]):
        origin = "bacteria"
        idx_list.append(0)
    if any([x in string for x in plants_check]):
        origin = "plant"
        idx_list.append(1)
    if any([x in string for x in fungi_check]):
        origin = "fungi"
        idx_list.append(2)
    if any([x in string for x in animal_check]):
        origin = "animal"
        idx_list.append(3)
    if any([x in string for x in marine_check]):
        origin = "marine"
        idx_list.append(4)
    
    if len(idx_list) == 0:
        origin = None
        idx = None
    elif len(idx_list) > 1:
        origin = "mix"
        idx = 5
    else:
        idx = idx_list[0]
        
    return origin, idx

#--------------------------------------------------------------------------#

def clean_smiles(
        smiles: str
        ) -> Chem.rdchem.Mol:
    """Uses the ChEMBL standardizer to parse a SMILES string
    
    Args:
        smiles:     SMILES of the molecule

    Returns:
        The correctly parsed RDKIT Mol object or None if structure cleaning
        failed
    """
    try:
        #built Mol object from SMILES
        mol = Chem.MolFromSmiles(smiles)

        #Pass mol through ChEMBL standardizer
        mol = standardizer.get_parent_mol(mol)[0]

        #Parse again back and forth through RDKIT
        #(sometimes ChEMBL messes up certain Mol properties)
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol),
                                 sanitize=True)
    except:
        mol = None

    return mol

#--------------------------------------------------------------------------#

def fix_dataframe(
        df: pd.DataFrame
        ) -> pd.DataFrame:
    """Processes raw COCONUT df for downstream analysis
    
    Args:
        df:     raw COCONUT df

    Returns:
        Processed dataframe with SMILES, organism, organism ID and
        Murko scaffold
    """


    output = df.copy()
    
    text = list(df["textTaxa"])
    smiles = list(df["sugar_free_smiles"])

    organism = [0]*len(df)
    class_id = [0]*len(df)
    mols = [0]*len(df)

    for i in range(len(df)):
        organism[i], class_id[i] = label_text(text[i])
        mols[i] = clean_smiles(smiles[i])

    output["organism"] = organism
    output["class_id"] = class_id
    output["mols"] = mols
    
    output = output.dropna()
    output = output[["sugar_free_smiles",
                     "organism",
                     "class_id",
                     "murko_framework",
                    ]]

    return output 





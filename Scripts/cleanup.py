import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
from chembl_structure_pipeline import standardizer


"""
Preprocessing routine adapted from:
https://github.com/reymond-group/Coconut-TMAP-SVM
"""

###########################################################################

def label_text(string):
    string = string.lower()
    
    bacteria_check = ["bacteria", "bacillus", "bacta", "bacterium"]
    plants_check = ["plants", "plant"]
    fungi_check = ["fungi", "fungus", "aspergillus"]
    animal_check = ["animal", "animals"]
    marine_check = ["marine"]

    if any([x in string for x in bacteria_check]):
        origin = "bacteria"
        idx = 0    
    elif any([x in string for x in plants_check]):
        origin = "plant"
        idx = 1
    elif any([x in string for x in fungi_check]):
        origin = "fungi"
        idx = 2
    elif any([x in string for x in animal_check]):
        origin = "animal"
        idx = 3
    elif any([x in string for x in marine_check]):
        origin = "marine"
        idx = 4
    else:
        origin = None
        idx = None
    
    return origin, idx
           

def clean_smiles(smiles):

    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = standardizer.get_parent_mol(mol)[0]
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol),
                                 sanitize=True)
    except:
        mol = None

    return mol


def fix_dataframe(df):
    
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
    output = output[["sugar_free_smiles", "organism",
                    "class_id"]]

    return output 

















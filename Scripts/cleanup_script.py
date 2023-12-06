"""COCONUT database cleaning script

Script to preprocess COCONUT into a .csv file that can be used for
further analysis. As long as the COCONUT_DB.sdf file is saved in the
../Data folder, it requires no user input.

Steps:
  1. Convert taxonomy information into known origin labels (i.e. "plant")
  2. Remove compounds that did not have appropriate taxonomy info
  3. Standardize SMILES and remove salts via ChEMBL pipeline standardizer
  4. Remove compounds that failed standardization
  5. Select only relevant columns (SMILES, origin, label, Murko scaffold)
  6. Analyse class statistics of cleaned dataset
  7. Save everything in ../Data
"""

from cleanup import *
from utils import *
from utils import get_statistics
from rdkit.Chem import PandasTools
import numpy as np
import argparse

###########################################################################

parser = argparse.ArgumentParser(description=__doc__,
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
args = parser.parse_args()

###########################################################################

def main():
  
  #initialize empty log list
  logs = []
  print("[cleanup]: Starting to load data")
  
  #load raw csv file
  coconut = PandasTools.LoadSDF("../Data/COCONUT_DB.sdf")
  logs.append("[cleanup]: Loaded " + str(len(coconut)) + " molecules") 
  print(logs[0])
  
  #remove compounds without taxonomy
  coconut = coconut.query("textTaxa != '[notax]'").copy()
  coconut = coconut.query("textTaxa != '[]'").copy()
  logs.append("[cleanup]: Found " + str(len(coconut)) + " molecules with taxonomy annotation")
  print(logs[1])
  
  #select relevant columns
  coconut = coconut[["sugar_free_smiles", "textTaxa", "murko_framework"]]
  
  #run taxonomy labelling and SMILES standardization
  coconut = fix_dataframe(coconut)
  logs.append("[cleanup]: Found " + str(len(coconut)) +
              " molecules with parseable taxonomy and correct structure")
  print(logs[2])
  
  #create numpy array with class ID information as one-hot encoding
  y = np.array(coconut["class_id"], dtype=np.int8)
  y = np.eye(6)[y]

  #get class statistics from the clean csv file
  stats = get_statistics(coconut)
  
  #print saving info and stats
  logs.append("[cleanup]: Saving data as ../Data/coconut.csv")
  logs.append("[cleanup]: Saving labels as ../Data/labels.pkl")
  logs.append("[cleanup]: Saving stats as ../Data/stats.csv")
  logs.append("[cleanup]: saving log as ../Data/log.txt")
  print(logs[3])
  print(logs[4])
  print(logs[5])
  print(logs[6])
  print(stats)

  #save everything
  coconut.to_csv("../Data/coconut.csv")
  stats.to_csv("../Data/stats.csv")
  with open("../Data/log.txt", 'w') as f:
    for line in logs:
            f.write(f"{line}\n")
  pickle_save(y, "../Data/labels.pkl", verbose=False) 
  



if __name__ == "__main__":             
  main()



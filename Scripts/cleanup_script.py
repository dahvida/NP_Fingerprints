from cleanup import *
from utils import get_statistics
from rdkit.Chem import PandasTools

###########################################################################

def main():

  logs = []
  print("[cleanup]: Starting to load data")
   
  coconut = PandasTools.LoadSDF("../Data/COCONUT_DB.sdf")
  logs.append("[cleanup]: Loaded " + str(len(coconut)) + " molecules") 
  print(logs[0])
  
  coconut = coconut.query("textTaxa != '[notax]'").copy()
  coconut = coconut.query("textTaxa != '[]'").copy()
  logs.append("[cleanup]: Found " + str(len(coconut)) + " molecules with taxonomy annotation")
  print(logs[1])
  
  coconut = coconut[["sugar_free_smiles", "textTaxa", "murko_framework"]]
  coconut = fix_dataframe(coconut)
  logs.append("[cleanup]: Found " + str(len(coconut)) +
              " molecules with parseable taxonomy and correct structure")
  print(logs[2])
  
  stats = get_statistics(coconut)
  
  logs.append("[cleanup]: Saving data as ../Data/clean_coconut.csv")
  logs.append("[cleanup]: Saving stats as ../Data/stats.csv")
  logs.append("[cleanup]: saving log as ../Data/log.txt")
  print(logs[3])
  print(logs[4])
  print(logs[5])
  print(stats)

  coconut.to_csv("../Data/clean_coconut.csv")
  stats.to_csv("../Data/stats.csv")
  with open("../Data/log.txt", 'w') as f:
    for line in logs:
            f.write(f"{line}\n")
  



if __name__ == "__main__":             
  main()



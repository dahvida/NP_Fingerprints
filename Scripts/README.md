# User notice
Unfortunately it was not possible to structure the package so that both jCompoundMapper and CDK fingerprints can be computed in the same session. This is due to a number of reasons:  
1. In order to use them, you have to use JPype to load them into Python since originally they were written in Java  
2. JPype does not allow starting separate JVM, nor restarting a JVM in the same session after closing one  
3. It is possible to load both of them into a single JVM, but then since a) both use CDK b) they require different versions of CDK c) both have stored CDK in exactly the same path of the Java class tree, conflicts will happen when trying to compute either set of fingerprints.  
In the future I plan to revise this, but for the moment either you import the fingerprints in `FP_calc.jmap_fps.py` or in `FP_calc.cdk_fps.py`. Any suggestion on how to overcome this issue, without reimplementing all fingerprints from scratch natively in Python, is welcome!  

# Folder structure

- [FP_calc:](clean_coconut.csv) Contains all functions and supporting files to compute
all fingerprints. Currently the package can calculate:  
   - MACCS  
   - Extended connectivity  
   - Functional class  
   - RDKIT  
   - Atom Pair  
   - Topological torsion  
   - Avalon  
   - MinHashed  
   - MAP4  
   - PUBCHEM   
   - Daylight   
   - ESTATE   
   - Klekota - Roth   
   - LINGO   
   - Depth First Search   
   - All Shortest Paths   
   - LSTAR   
   - RAD2D   
   - Pharmacophoric pairs   
   - Pharmacophoric triplets  

- [cleanup.py:](cleanup.py) Contains utility functions used in `cleanup_script`.  

- [cleanup_script.py:](cleanup_script.py) Preprocesses the raw COCONUT database. For further
information, open the file or use the `--help` flag from the command line.  

- [clf.py:](clf.py) Contains utility functions used in `clf_script`.  

- [clf_script.py:](clf_script.py) Calculates classification metrics for all classes
using all precomputed fingerprints. For further information, open the file or use the `--help` flag from the command line.  

- [fp_script.py:](fp_script.py) Computes the chosen fingerprint set. For further information, open the file or use the `--help` flag from the command line.  

- [sim_search.py:](sim_search.py) Contains utility functions used in `sim_search_script`.  

- [sim_search_script.py:](sim_search_script.py) Runs the similarity search using all precomputed fingerprints for the COCONUT dataset. For further information,
 open the file or use the `--help` flag from the command line.  

- [drug_rep_hub.py:](drug_rep_hub.py) Runs the similarity search and correlation matrix analysis for the Drug Repurposing Hub dataset.

- [utils.py:](utils.py) Contains utility functions used throughout the repository.  

